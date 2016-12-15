      MODULE F77_LAPACK
!
!  -- LAPACK95 interface driver routine (version 3.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     September, 2000
!
 
      INTERFACE LA_LANGB

       FUNCTION DLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLANGB
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDAB, N, KL, KU
         REAL(WP), INTENT(IN) :: AB( LDAB, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION DLANGB

       FUNCTION ZLANGB( NORM, N, KL, KU, AB, LDAB, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: ZLANGB
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDAB, N, KL, KU
         COMPLEX(WP), INTENT(IN) :: AB( LDAB, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION ZLANGB

       END INTERFACE
       
       INTERFACE LA_TGSEN

      SUBROUTINE DTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, &
     &                   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL,   &
     &                   PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      LOGICAL, INTENT(IN) :: WANTQ, WANTZ
      INTEGER, INTENT(IN) :: IJOB, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, N
      INTEGER, INTENT(OUT) :: INFO, M, IWORK(LIWORK)
      REAL(WP), INTENT(OUT) :: PL, PR
      LOGICAL, INTENT(IN) :: SELECT(*)
      REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*), Z(LDZ,*)
      REAL(WP), INTENT(OUT) :: ALPHAI(*), ALPHAR(*), BETA(*), DIF(2),   &
     &                         WORK(LWORK)
      END SUBROUTINE DTGSEN


      SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB, &
     &                   ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,   &
     &                   WORK, LWORK, IWORK, LIWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      LOGICAL, INTENT(IN) :: WANTQ, WANTZ
      INTEGER, INTENT(IN) :: IJOB, LDA, LDB, LDQ, LDZ, LIWORK, LWORK, N
      INTEGER, INTENT(OUT) :: INFO, M, IWORK(LIWORK)
      REAL(WP), INTENT(OUT) :: PL, PR
      LOGICAL, INTENT(IN) :: SELECT(*)
      REAL(WP), INTENT(OUT) :: DIF(2)
      COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),       &
     &                              Z(LDZ,*)
      COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), WORK(LWORK)
      END SUBROUTINE ZTGSEN

       END INTERFACE

       
       INTERFACE LA_TGSNA

      SUBROUTINE DTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,    &
     &                   LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,    &
     &                   IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, JOB
      INTEGER, INTENT(IN) :: LDA, LDB, LDVL, LDVR, LWORK, MM, N
      INTEGER, INTENT(OUT) :: INFO, M, IWORK(*)
      LOGICAL, INTENT(IN) :: SELECT(*)
      REAL(WP), INTENT(OUT) :: DIF(*), S(*)
      REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*), VL(LDVL,*),           &
     &                        VR(LDVR,*)
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DTGSNA

      SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,    &
     &                   LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,    &
     &                   IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, JOB
      INTEGER, INTENT(IN) :: LDA, LDB, LDVL, LDVR, LWORK, MM, N
      INTEGER, INTENT(OUT) :: INFO, M, IWORK(*)
      LOGICAL, INTENT(IN) :: SELECT(*)
      REAL(WP), INTENT(OUT) :: DIF(*), S(*)
      COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*), VL(LDVL,*),        &
     &                           VR(LDVR,*)
      COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZTGSNA

       END INTERFACE
       
       INTERFACE LA_TGSYL

      SUBROUTINE DTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,  &
     &                   LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,  &
     &                   IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: TRANS
      INTEGER, INTENT(IN) :: IJOB, LDA, LDB, LDC, LDD, LDE, LDF, LWORK, &
     &                       M, N
      INTEGER, INTENT(OUT) :: INFO, IWORK(*)
      REAL(WP), INTENT(OUT) :: DIF, SCALE
      REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*), D(LDD,*), E(LDF,*)
      REAL(WP), INTENT(INOUT) :: C(LDC,*), F(LDF,*)
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DTGSYL

      SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,  &
     &                   LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,  &
     &                   IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: TRANS
      INTEGER, INTENT(IN) :: IJOB, LDA, LDB, LDC, LDD, LDE, LDF, LWORK, &
     &                       M, N
      INTEGER, INTENT(OUT) :: INFO, IWORK(*)
      REAL(WP), INTENT(OUT) :: DIF, SCALE
      COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*), D(LDD,*), E(LDF,*)
      COMPLEX(WP), INTENT(INOUT) :: C(LDC,*), F(LDF,*)
      COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZTGSYL

       END INTERFACE
       
       INTERFACE LA_TGEXC

         SUBROUTINE DTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, &
     &                      LDZ, IFST, ILST, WORK, LWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      LOGICAL, INTENT(IN) :: WANTQ, WANTZ
      INTEGER, INTENT(IN) :: LDA, LDB, LDQ, LDZ, LWORK, N
      INTEGER, INTENT(INOUT) :: IFST, ILST
      INTEGER, INTENT(OUT) :: INFO
      REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*), Z(LDZ,*)
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DTGEXC


         SUBROUTINE ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, &
     &                      LDZ, IFST, ILST, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      LOGICAL, INTENT(IN) :: WANTQ, WANTZ
      INTEGER, INTENT(IN) ::  LDA, LDB, LDQ, LDZ, N
      INTEGER, INTENT(INOUT) :: IFST, ILST
      INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),    &
     &                                 Z(LDZ,*)
      END SUBROUTINE ZTGEXC

       END INTERFACE

       INTERFACE LA_BDSDC

         SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q,  &
     &                      IQ, WORK, IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: COMPQ, UPLO
      INTEGER, INTENT(IN) :: LDU, LDVT, N
      INTEGER, INTENT(OUT) :: INFO, IQ( * ), IWORK( * )
      REAL(WP), INTENT(INOUT) :: D( * ), E( * )
      REAL(WP), INTENT(OUT) :: Q(*), U(LDU,*), VT(LDVT,*), WORK(*)
      END SUBROUTINE DBDSDC

       END INTERFACE

       INTERFACE LA_STEGR

         SUBROUTINE DSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,       &
     &                      ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,  &
     &                      IWORK, LIWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP  
      CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE
      INTEGER, INTENT(IN) :: IL, IU, LDZ, LIWORK, LWORK, N
      INTEGER, INTENT(OUT) :: INFO, M
      INTEGER, INTENT(OUT) :: ISUPPZ( * ), IWORK(LIWORK)
      REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
      REAL(WP), INTENT(INOUT) :: D( * ), E( * )
      REAL(WP), INTENT(IN) :: W( * )
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      REAL(WP), INTENT(OUT) :: Z( LDZ, * )
      END SUBROUTINE DSTEGR
        
         SUBROUTINE ZSTEGR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU,       &
     &                      ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,  &
     &                      IWORK, LIWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP  
      CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE
      INTEGER, INTENT(IN) :: IL, IU, LDZ, LIWORK, LWORK, N
      INTEGER, INTENT(OUT) :: INFO, M
      INTEGER, INTENT(OUT) :: ISUPPZ( * ), IWORK(LIWORK)
      REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
      REAL(WP), INTENT(INOUT) :: D( * ), E( * )
      REAL(WP), INTENT(IN) :: W( * )
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      COMPLEX(WP), INTENT(OUT) :: Z( LDZ, * )
      END SUBROUTINE ZSTEGR
        
       END INTERFACE

       INTERFACE LA_ORMRZ

         SUBROUTINE DORMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C,    &
     &                      LDC, WORK, LWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
      INTEGER, INTENT(IN) :: K, L, LDA, LDC, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      REAL(WP), INTENT(IN) :: A( LDA, * ), TAU( * )
      REAL(WP), INTENT(INOUT) :: C( LDC, * )
      REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMRZ

       END INTERFACE


       INTERFACE LA_UNMRZ

         SUBROUTINE ZUNMRZ( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C,    &
     &                      LDC, WORK, LWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
      INTEGER, INTENT(IN) :: K, L, LDA, LDC, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      COMPLEX(WP), INTENT(IN) :: A( LDA, * ), TAU( * )
      COMPLEX(WP), INTENT(INOUT) :: C( LDC, * )
      COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMRZ

       END INTERFACE

       INTERFACE LA_TZRZF

         SUBROUTINE DTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      INTEGER, INTENT(IN) :: LDA, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      REAL(WP), INTENT(INOUT) :: A( LDA, * )
      REAL(WP), INTENT(OUT) :: TAU( * ), WORK(LWORK)
      END SUBROUTINE DTZRZF

         SUBROUTINE ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      INTEGER, INTENT(IN) :: LDA, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      COMPLEX(WP), INTENT(INOUT) :: A( LDA, * )
      COMPLEX(WP), INTENT(OUT) :: TAU( * ), WORK(LWORK)
      END SUBROUTINE ZTZRZF

       END INTERFACE

       INTERFACE LA_GEQP3

         SUBROUTINE DGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK,       &
     &                      INFO )
      USE LA_PRECISION, ONLY: WP => DP
      INTEGER, INTENT(IN) :: LDA, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(INOUT) :: JPVT( * )
      REAL(WP), INTENT(INOUT) :: A( LDA, * )
      REAL(WP), INTENT(OUT) :: TAU( * ), WORK(LWORK)
      END SUBROUTINE DGEQP3


         SUBROUTINE ZGEQP3( M, N, A, LDA, JPVT, TAU, WORK, LWORK, RWORK,&
     &                      INFO )
      USE LA_PRECISION, ONLY: WP => DP
      INTEGER, INTENT(IN) :: LDA, LWORK, M, N
      INTEGER, INTENT(OUT) :: INFO
      INTEGER, INTENT(INOUT) :: JPVT( * )
      COMPLEX(WP), INTENT(INOUT) :: A( LDA, * )
      COMPLEX(WP), INTENT(OUT) :: TAU( * ), WORK(LWORK)
      REAL(WP), INTENT(OUT) ::  RWORK( * )
      END SUBROUTINE ZGEQP3

       END INTERFACE

       INTERFACE LA_GESDD


         SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,    &
     &                      WORK, LWORK, IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: JOBZ
      INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
      INTEGER, INTENT(OUT) :: INFO
      REAL(WP), INTENT(OUT) :: S(*)
      REAL(WP), INTENT(INOUT) :: A(LDA,*)
      REAL(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
      INTEGER :: IWORK(*)
      END SUBROUTINE DGESDD


        SUBROUTINE ZGESDD(  JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT,    &
     &                      WORK, LWORK, RWORK, IWORK, INFO )
      USE LA_PRECISION, ONLY: WP => DP
      CHARACTER(LEN=1), INTENT(IN) :: JOBZ
      INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
      INTEGER, INTENT(OUT) :: INFO
      REAL(WP), INTENT(OUT) :: S(*)
      REAL(WP) :: RWORK(*)
      COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
      COMPLEX(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
      INTEGER :: IWORK(*)
      END SUBROUTINE  ZGESDD
      END INTERFACE       


      INTERFACE LA_GGRQF

      SUBROUTINE DGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, LWORK, M, N, P
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: TAUA(*), TAUB(*), WORK(*)
      END SUBROUTINE DGGRQF

      SUBROUTINE ZGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, LWORK, M, N, P
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: TAUA(*), TAUB(*), WORK(*)
      END SUBROUTINE ZGGRQF

      END INTERFACE

      INTERFACE LA_GGQRF

      SUBROUTINE DGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, LWORK, M, N, P
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: TAUA(*), TAUB(*), WORK(*)
      END SUBROUTINE DGGQRF

      SUBROUTINE ZGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, LWORK, M, N, P
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: TAUA(*), TAUB(*), WORK(*)
      END SUBROUTINE ZGGQRF

      END INTERFACE

      INTERFACE LA_DISNA

      SUBROUTINE DDISNA( JOB, M, N, D, SEP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB
         INTEGER, INTENT(IN) :: M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(OUT) :: SEP(*)
      END SUBROUTINE DDISNA

      END INTERFACE

      INTERFACE LA_TGSJA

      SUBROUTINE DTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,    &
     &                   LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,  &
     &                   Q, LDQ, WORK, NCYCLE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBQ, JOBU, JOBV
         INTEGER, INTENT(IN) :: K, L, LDA, LDB, LDQ, LDU, LDV, M, N,    &
     &                          NCYCLE, P
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TOLA, TOLB
         REAL(WP), INTENT(OUT) :: ALPHA(*), BETA(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),       &
     &                              U(LDU,*), V(LDV,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTGSJA

      SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,    &
     &                   LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,  &
     &                   Q, LDQ, WORK, NCYCLE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBQ, JOBU, JOBV
         INTEGER, INTENT(IN) :: K, L, LDA, LDB, LDQ, LDU, LDV, M, N,    &
     &                          NCYCLE, P
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TOLA, TOLB
         REAL(WP), INTENT(OUT) :: ALPHA(*), BETA(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),    &
     &                                 U(LDU,*), V(LDV,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTGSJA

      END INTERFACE

      INTERFACE LA_GGSVP

      SUBROUTINE DGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,     &
     &                   TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,      &
     &                   IWORK, TAU, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBQ, JOBU, JOBV
         INTEGER, INTENT(IN) :: LDA, LDB, LDQ, LDU, LDV, M, N, P
         INTEGER, INTENT(OUT) :: INFO, K, L, IWORK(*)
         REAL(WP), INTENT(IN) :: TOLA, TOLB
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: Q(LDQ,*), TAU(*), U(LDU,*), V(LDV,*), &
     &                            WORK(*)
      END SUBROUTINE DGGSVP

      SUBROUTINE ZGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,     &
     &                   TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,      &
     &                   IWORK, RWORK, TAU, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBQ, JOBU, JOBV
         INTEGER, INTENT(IN) :: LDA, LDB, LDQ, LDU, LDV, M, N, P
         INTEGER, INTENT(OUT) :: INFO, K, L, IWORK(*)
         REAL(WP), INTENT(IN) :: TOLA, TOLB
         REAL(WP), INTENT(IN) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: Q(LDQ,*), TAU(*), U(LDU,*),        &
     &                               V(LDV,*), WORK(*)
      END SUBROUTINE ZGGSVP

      END INTERFACE

      INTERFACE LA_TGEVC

      SUBROUTINE DTGEVC( SIDE, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,   &
     &                   LDVL, VR, LDVR, MM, M, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, SIDE
         INTEGER, INTENT(IN) :: LDA, LDB, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: VL(LDVL,*), VR(LDVR,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTGEVC

      SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,   &
     &                   LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, SIDE
         INTEGER, INTENT(IN) :: LDA, LDB, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: VL(LDVL,*), VR(LDVR,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTGEVC

      END INTERFACE

      INTERFACE LA_HGEQZ

      SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,&
     &                   ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,    &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, COMPZ, JOB
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDB, LDQ, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),       &
     &                              Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*), WORK(LWORK)
      END SUBROUTINE DHGEQZ

      SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,&
     &                   ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,      &
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, COMPZ, JOB
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDB, LDQ, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RWORK( * )
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),    &
     &                                 Z(LDZ,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), WORK(LWORK)
      END SUBROUTINE ZHGEQZ

      END INTERFACE

      INTERFACE LA_GGBAK

      SUBROUTINE DGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,  &
     &                   LDV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB, SIDE
         INTEGER, INTENT(IN) :: IHI, ILO, LDV, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: LSCALE(*), RSCALE(*)
         REAL(WP), INTENT(INOUT) :: V(LDV,*)
      END SUBROUTINE DGGBAK

      SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,  &
     &                   LDV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB, SIDE
         INTEGER, INTENT(IN) :: IHI, ILO, LDV, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: LSCALE(*), RSCALE(*)
         COMPLEX(WP), INTENT(INOUT) :: V(LDV,*)
      END SUBROUTINE ZGGBAK

      END INTERFACE

      INTERFACE LA_GGBAL

      SUBROUTINE DGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,      &
     &                   RSCALE, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB
         INTEGER, INTENT(IN) :: LDA, LDB, N
         INTEGER, INTENT(OUT) :: IHI, ILO, INFO
         REAL(WP), INTENT(OUT) :: LSCALE(*), RSCALE(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DGGBAL

      SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,      &
     &                   RSCALE, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB
         INTEGER, INTENT(IN) :: LDA, LDB, N
         INTEGER, INTENT(OUT) :: IHI, ILO, INFO
         REAL(WP), INTENT(OUT) :: LSCALE(*), RSCALE(*), WORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE ZGGBAL

      END INTERFACE

      INTERFACE LA_GGHRD

      SUBROUTINE DGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,  &
     &                   LDQ, Z, LDZ, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, COMPZ
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDB, LDQ, LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),       &
     &                              Z(LDZ,*)
      END SUBROUTINE DGGHRD

      SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,  &
     &                   LDQ, Z, LDZ, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, COMPZ
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDB, LDQ, LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), Q(LDQ,*),    &
     &                                 Z(LDZ,*)
      END SUBROUTINE ZGGHRD

      END INTERFACE

      INTERFACE LA_PBSTF

      SUBROUTINE DPBSTF( UPLO, N, KD, AB, LDAB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) ::UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE DPBSTF

      SUBROUTINE ZPBSTF( UPLO, N, KD, AB, LDAB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) ::UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB( LDAB, * )
      END SUBROUTINE ZPBSTF

      END INTERFACE

      INTERFACE LA_SBGST

      SUBROUTINE DSBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,  &
     &                   LDX, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT
         INTEGER, INTENT(IN) :: KA, KB, LDAB, LDBB, LDX, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: BB(LDBB,*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: WORK(*), X(LDX,*)
      END SUBROUTINE DSBGST

      END INTERFACE

      INTERFACE LA_HBGST

      SUBROUTINE ZHBGST( VECT, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, X,  &
     &                   LDX, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT
         INTEGER, INTENT(IN) :: KA, KB, LDAB, LDBB, LDX, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(IN) :: BB(LDBB,*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), X(LDX,*)
      END SUBROUTINE ZHBGST

      END INTERFACE

      INTERFACE LA_SPGST

      SUBROUTINE DSPGST( ITYPE, UPLO, N, AP, BP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: ITYPE, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: BP(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE DSPGST

      END INTERFACE

      INTERFACE LA_HPGST

      SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: ITYPE, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: BP(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE ZHPGST

      END INTERFACE

      INTERFACE LA_BDSQR

      SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,    &
     &                   LDU, C, LDC, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
         INTEGER, INTENT(OUT) :: INFO
         REAL, INTENT(INOUT) :: D(*), E(*)
         REAL, INTENT(OUT) :: RWORK(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*), U(LDU,*), VT(LDVT,*)
      END SUBROUTINE DBDSQR

      SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,    &
     &                   LDU, C, LDC, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDC, LDU, LDVT, N, NCC, NCVT, NRU
         INTEGER, INTENT(OUT) :: INFO
         REAL, INTENT(INOUT) :: D(*), E(*)
         REAL, INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*), U(LDU,*), VT(LDVT,*)
      END SUBROUTINE ZBDSQR

      END INTERFACE

      INTERFACE LA_ORMBR

      SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,    &
     &                   LDC, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, VECT
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMBR

      END INTERFACE

      INTERFACE LA_UNMBR

      SUBROUTINE ZUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,    &
     &                   LDC, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, VECT
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMBR

      END INTERFACE

      INTERFACE LA_ORGBR

      SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK,       &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: VECT
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGBR

      END INTERFACE

      INTERFACE LA_UNGBR

      SUBROUTINE ZUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK,       &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: VECT
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGBR

      END INTERFACE

      INTERFACE LA_GBBRD

      SUBROUTINE DGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q,    &
     &                   LDQ, PT, LDPT, C, LDC, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: VECT
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC
         INTEGER, INTENT(OUT) :: INFO
         REAL, INTENT(OUT) :: D(*), E(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), C(LDC,*)
         REAL(WP), INTENT(OUT) :: PT(LDPT,*), Q(LDQ,*), WORK(*)
      END SUBROUTINE DGBBRD

      SUBROUTINE ZGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q,    &
     &                   LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: VECT
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC
         INTEGER, INTENT(OUT) :: INFO
         REAL, INTENT(OUT) :: D(*), E(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: PT(LDPT,*), Q(LDQ,*), WORK(*)
      END SUBROUTINE ZGBBRD

      END INTERFACE

      INTERFACE LA_GEBRD

      SUBROUTINE DGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,   &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAUP(*), TAUQ(*), WORK(LWORK)
      END SUBROUTINE DGEBRD

      SUBROUTINE ZGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,   &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAUP(*), TAUQ(*), WORK(LWORK)
      END SUBROUTINE ZGEBRD

      END INTERFACE

      INTERFACE LA_TRSEN

      SUBROUTINE DTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI, &
     &                   M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, JOB
         INTEGER, INTENT(IN) :: LDQ, LDT, LWORK, N, LIWORK
         INTEGER, INTENT(OUT) :: INFO, M, IWORK(LIWORK)
         REAL(WP), INTENT(OUT) :: S, SEP
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(INOUT) :: Q(LDQ,*), T(LDT,*)
         REAL(WP), INTENT(IN) :: WR(*), WI(*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DTRSEN

      SUBROUTINE ZTRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, W, M, S,&
     &                   SEP, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ, JOB
         INTEGER, INTENT(IN) :: LDQ, LDT, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, M
         REAL(WP), INTENT(OUT) :: S, SEP
         LOGICAL, INTENT(IN) :: SELECT(*)
         COMPLEX(WP), INTENT(INOUT) :: Q(LDQ,*), T(LDT,*)
         COMPLEX(WP), INTENT(IN) :: W(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZTRSEN

      END INTERFACE

      INTERFACE LA_TRSNA

      SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,  &
     &                   LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,      &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, JOB
         INTEGER, INTENT(IN) :: LDT, LDVL, LDVR, LDWORK, MM, N
         INTEGER, INTENT(OUT) :: INFO, M, IWORK(*)
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(OUT) :: S(*), SEP(*)
         REAL(WP), INTENT(IN) :: T(LDT,*), VL(LDVL,*), VR(LDVR,*)
         REAL(WP), INTENT(OUT) :: WORK(LDWORK,*)
      END SUBROUTINE DTRSNA

      SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,  &
     &                   LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK,      &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, JOB
         INTEGER, INTENT(IN) :: LDT, LDVL, LDVR, LDWORK, MM, N
         INTEGER, INTENT(OUT) :: INFO, M
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(OUT) :: RWORK(*), S(*), SEP(*)
         COMPLEX(WP), INTENT(IN) :: T(LDT,*), VL(LDVL,*), VR(LDVR,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LDWORK,*)
      END SUBROUTINE ZTRSNA

      END INTERFACE

      INTERFACE LA_TRSYL

      SUBROUTINE DTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,   &
     &                   LDC, SCALE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANA, TRANB
         INTEGER, INTENT(IN) :: ISGN, LDA, LDB, LDC, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: SCALE
         REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
      END SUBROUTINE DTRSYL

      SUBROUTINE ZTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,   &
     &                   LDC, SCALE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANA, TRANB
         INTEGER, INTENT(IN) :: ISGN, LDA, LDB, LDC, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: SCALE
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
      END SUBROUTINE ZTRSYL

      END INTERFACE

      INTERFACE LA_TREXC

      SUBROUTINE DTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK,    &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ
         INTEGER, INTENT(IN) :: IFST, ILST, LDQ, LDT, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: Q(LDQ,*), T(LDT,*), WORK(*)
      END SUBROUTINE DTREXC

      SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPQ
         INTEGER, INTENT(IN) :: IFST, ILST, LDQ, LDT, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: Q(LDQ,*), T(LDT,*)
      END SUBROUTINE ZTREXC

      END INTERFACE

      INTERFACE LA_TREVC

      SUBROUTINE DTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
     &                   LDVR, MM, M, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, SIDE
         INTEGER, INTENT(IN) :: LDT, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M
         LOGICAL, INTENT(INOUT) :: SELECT(*)
         REAL(WP), INTENT(IN) :: T(LDT,*)
         REAL(WP), INTENT(INOUT) :: VL(LDVL,*), VR(LDVR,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTREVC

      SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
     &                   LDVR, MM, M, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: HOWMNY, SIDE
         INTEGER, INTENT(IN) :: LDT, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M
         LOGICAL, INTENT(INOUT) :: SELECT(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: T(LDT,*), VL(LDVL,*), VR(LDVR,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTREVC

      END INTERFACE

      INTERFACE LA_HSEIN

      SUBROUTINE DHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, WR, WI,&
     &                   VL, LDVL, VR, LDVR, MM, M, WORK, IFAILL,       &
     &                   IFAILR, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: EIGSRC, INITV, SIDE
         INTEGER, INTENT(IN) :: LDH, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M, IFAILL(*), IFAILR(*)
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(INOUT) :: WR(*), WI(*)
         REAL(WP), INTENT(IN) :: H(LDH,*)
         REAL(WP), INTENT(INOUT) :: VL(LDVL,*), VR(LDVR,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DHSEIN

      SUBROUTINE ZHSEIN( SIDE, EIGSRC, INITV, SELECT, N, H, LDH, W, VL, &
     &                   LDVL, VR, LDVR, MM, M, WORK, RWORK, IFAILL,    &
     &                   IFAILR, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: EIGSRC, INITV, SIDE
         INTEGER, INTENT(IN) :: LDH, LDVL, LDVR, MM, N
         INTEGER, INTENT(OUT) :: INFO, M, IFAILL(*), IFAILR(*)
         LOGICAL, INTENT(IN) :: SELECT(*)
         REAL(WP), INTENT(OUT) :: RWORK( * )
         COMPLEX(WP), INTENT(IN) :: H(LDH,*)
         COMPLEX(WP), INTENT(INOUT) :: VL(LDVL,*), VR(LDVR,*), W(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHSEIN

      END INTERFACE

      INTERFACE LA_HSEQR

      SUBROUTINE DHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, WR, WI, Z,    &
     &                   LDZ, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ, JOB
         INTEGER, INTENT(IN) :: IHI, ILO, LDH, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: WR(*), WI(*)
         REAL(WP), INTENT(INOUT) :: H(LDH,*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DHSEQR

      SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,    &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ, JOB
         INTEGER, INTENT(IN) :: IHI, ILO, LDH, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: H(LDH,*), Z(LDZ,*)
         COMPLEX(WP), INTENT(OUT) :: W(*), WORK(LWORK)
      END SUBROUTINE ZHSEQR

      END INTERFACE

      INTERFACE LA_ORMHR

      SUBROUTINE DORMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,   &
     &                   LDC, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1),  INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMHR

      END INTERFACE

      INTERFACE LA_UNMHR

      SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,   &
     &                   LDC, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1),  INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMHR

      END INTERFACE

      INTERFACE LA_ORGHR

      SUBROUTINE DORGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGHR

      END INTERFACE

      INTERFACE LA_UNGHR

      SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGHR

      END INTERFACE

      INTERFACE LA_GEBAK

      SUBROUTINE DGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,      &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB, SIDE
         INTEGER, INTENT(IN) :: IHI, ILO, LDV, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: SCALE(*)
         REAL(WP), INTENT(INOUT) :: V(LDV,*)
      END SUBROUTINE DGEBAK

      SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,      &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB, SIDE
         INTEGER, INTENT(IN) :: IHI, ILO, LDV, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: SCALE(*)
         COMPLEX(WP), INTENT(INOUT) :: V(LDV,*)
      END SUBROUTINE ZGEBAK

      END INTERFACE

      INTERFACE LA_GEBAL

      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: IHI, ILO, INFO
         REAL(WP), INTENT(OUT) :: SCALE(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE DGEBAL

      SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOB
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: IHI, ILO, INFO
         REAL(WP), INTENT(OUT) :: SCALE(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE ZGEBAL

      END INTERFACE

      INTERFACE LA_GEHRD

      SUBROUTINE DGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DGEHRD

      SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: IHI, ILO, LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZGEHRD

      END INTERFACE

      INTERFACE LA_PTEQR

      SUBROUTINE DPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: INFO, LDZ, N
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(INOUT) :: Z(LDZ,*)
      END SUBROUTINE DPTEQR

      SUBROUTINE ZPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: INFO, LDZ, N
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(INOUT) :: Z(LDZ,*)
      END SUBROUTINE ZPTEQR

      END INTERFACE

      INTERFACE LA_STEIN

      SUBROUTINE DSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,   &
     &                   IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDZ, M, N, IBLOCK(*), ISPLIT(*)
         INTEGER, INTENT(OUT) :: INFO, IFAIL(*), IWORK(*)
         REAL(WP), INTENT(IN) :: D(*), E(*), W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(OUT) :: Z( LDZ, * )
      END SUBROUTINE DSTEIN

      SUBROUTINE ZSTEIN( N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,   &
     &                   IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDZ, M, N, IBLOCK(*), ISPLIT(*)
         INTEGER, INTENT(OUT) :: INFO, IFAIL(*), IWORK(*)
         REAL(WP), INTENT(IN) :: D(*), E(*), W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z( LDZ, * )
      END SUBROUTINE ZSTEIN

      END INTERFACE

      INTERFACE LA_STEBZ

      SUBROUTINE DSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
     &                   M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,     &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: ORDER, RANGE
         INTEGER, INTENT(IN) :: IL, IU, M, N
         INTEGER, INTENT(OUT) :: INFO, NSPLIT, IBLOCK(*), ISPLIT(*),    &
     &                           IWORK(*) 
         REAL(WP), INTENT(IN) :: ABSTOL, VL, VU, D(*), E(*)
         REAL(WP), INTENT(OUT) :: W(*), WORK(*)
      END SUBROUTINE DSTEBZ

      END INTERFACE

      INTERFACE LA_STEDC

      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,    &
     &                   LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: LDZ, LIWORK, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(LIWORK)
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(INOUT) :: Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DSTEDC

      SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK,    &
     &                   LRWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: LDZ, LIWORK, LRWORK, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(LIWORK)
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: RWORK(LRWORK)
         COMPLEX(WP), INTENT(INOUT) :: Z(LDZ,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZSTEDC

      END INTERFACE

      INTERFACE LA_STERF

      SUBROUTINE DSTERF( N, D, E, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
      END SUBROUTINE DSTERF

      END INTERFACE

      INTERFACE LA_STEQR

      SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(INOUT) :: Z(LDZ,*)
      END SUBROUTINE DSTEQR

      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: COMPZ
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(INOUT) :: Z(LDZ,*)
      END SUBROUTINE ZSTEQR

      END INTERFACE

      INTERFACE LA_OPMTR

      SUBROUTINE DOPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDC, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DOPMTR

      END INTERFACE

      INTERFACE LA_UPMTR

      SUBROUTINE ZUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDC, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZUPMTR

      END INTERFACE

      INTERFACE LA_OPGTR

      SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDQ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*), TAU(*)
         REAL(WP), INTENT(OUT) :: Q(LDQ,*), WORK(*)
      END SUBROUTINE DOPGTR

      END INTERFACE

      INTERFACE LA_UPGTR

      SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDQ, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*), TAU(*)
         COMPLEX(WP), INTENT(OUT) :: Q(LDQ,*), WORK(*)
      END SUBROUTINE ZUPGTR

      END INTERFACE

      INTERFACE LA_ORMTR

      SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,  &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
      END SUBROUTINE DORMTR

      END INTERFACE

      INTERFACE LA_UNMTR

      SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,  &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
      END SUBROUTINE ZUNMTR

      END INTERFACE

      INTERFACE LA_SBTRD

      SUBROUTINE DSBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,     &
     &                   WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT
         INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), Q(LDQ,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSBTRD

      END INTERFACE

      INTERFACE LA_HBTRD

      SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,     &
     &                   WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, VECT
         INTEGER, INTENT(IN) :: KD, LDAB, LDQ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), Q(LDQ,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHBTRD

      END INTERFACE

      INTERFACE LA_SPTRD

      SUBROUTINE DSPTRD( UPLO, N, AP, D, E, TAU, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: TAU(*)
      END SUBROUTINE DSPTRD

      END INTERFACE

      INTERFACE LA_HPTRD

      SUBROUTINE ZHPTRD( UPLO, N, AP, D, E, TAU, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*)
      END SUBROUTINE ZHPTRD

      END INTERFACE

      INTERFACE LA_TZRQF

      SUBROUTINE DTZRQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DTZRQF

      SUBROUTINE ZTZRQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZTZRQF

      END INTERFACE

      INTERFACE LA_ORMRQ

      SUBROUTINE DORMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMRQ

      END INTERFACE

      INTERFACE LA_UNMRQ

      SUBROUTINE ZUNMRQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMRQ

      END INTERFACE

      INTERFACE LA_ORGRQ

      SUBROUTINE DORGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGRQ

      END INTERFACE

      INTERFACE LA_UNGRQ

      SUBROUTINE ZUNGRQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGRQ

      END INTERFACE

      INTERFACE LA_GERQF

      SUBROUTINE DGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DGERQF

      SUBROUTINE ZGERQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZGERQF

      END INTERFACE

      INTERFACE LA_ORMQL

      SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMQL

      END INTERFACE

      INTERFACE LA_UNMQL

      SUBROUTINE ZUNMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMQL

      END INTERFACE

      INTERFACE LA_ORGQL

      SUBROUTINE DORGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGQL

      END INTERFACE

      INTERFACE LA_UNGQL

      SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGQL

      END INTERFACE

      INTERFACE LA_GEQLF

      SUBROUTINE DGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DGEQLF

      SUBROUTINE ZGEQLF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZGEQLF

      END INTERFACE

      INTERFACE LA_ORMLQ

      SUBROUTINE DORMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMLQ

      END INTERFACE

      INTERFACE LA_UNMLQ

      SUBROUTINE ZUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMLQ

      END INTERFACE

      INTERFACE LA_ORGLQ

      SUBROUTINE DORGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGLQ

      END INTERFACE

      INTERFACE LA_UNGLQ

      SUBROUTINE ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGLQ

      END INTERFACE

      INTERFACE LA_GELQF

      SUBROUTINE DGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DGELQF

      SUBROUTINE ZGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZGELQF

      END INTERFACE

      INTERFACE LA_ORMQR

      SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         REAL(WP), INTENT(INOUT) :: C(LDC,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORMQR

      END INTERFACE

      INTERFACE LA_UNMQR

      SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,     &
     &                   WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: SIDE, TRANS
         INTEGER, INTENT(IN) :: K, LDA, LDC, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: C(LDC,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNMQR

      END INTERFACE

      INTERFACE LA_ORGQR

      SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGQR

      END INTERFACE

      INTERFACE LA_UNGQR

      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: K, LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGQR

      END INTERFACE

      INTERFACE LA_GEQRF

      SUBROUTINE DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DGEQRF

      SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, M, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZGEQRF

      END INTERFACE

      INTERFACE LA_GEQPF

      SUBROUTINE DGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(*)
      END SUBROUTINE DGEQPF

      SUBROUTINE ZGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(*)
      END SUBROUTINE ZGEQPF

      END INTERFACE

      INTERFACE LA_TBRFS

      SUBROUTINE DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,   &
     &                   LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: AB(LDAB,*), B(LDB,*), X(LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTBRFS

      SUBROUTINE ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,   &
     &                   LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB(LDAB,*), B(LDB,*), X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTBRFS

      MODULE PROCEDURE DTBRFS1
      MODULE PROCEDURE ZTBRFS1

      END INTERFACE

      INTERFACE LA_TBCON

      SUBROUTINE DTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,&
     &                   IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(IN) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTBCON

      SUBROUTINE ZTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,&
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB(LDAB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTBCON

      END INTERFACE

      INTERFACE LA_TBTRS

      SUBROUTINE DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,   &
     &                   LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AB(LDAB,*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DTBTRS

      SUBROUTINE ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,   &
     &                   LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AB(LDAB,*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZTBTRS

      MODULE PROCEDURE DTBTRS1
      MODULE PROCEDURE ZTBTRS1

      END INTERFACE

      INTERFACE LA_TPTRI

      SUBROUTINE DTPTRI( UPLO, DIAG, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP( * )
      END SUBROUTINE DTPTRI

      SUBROUTINE ZTPTRI( UPLO, DIAG, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP( * )
      END SUBROUTINE ZTPTRI

      END INTERFACE

      INTERFACE LA_TPRFS

      SUBROUTINE DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,&
     &                   FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: AP(*), B(LDB,*), X(LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTPRFS

      SUBROUTINE ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,&
     &                   FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AP(*), B(LDB,*), X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTPRFS

      MODULE PROCEDURE DTPRFS1
      MODULE PROCEDURE ZTPRFS1

      END INTERFACE

      INTERFACE LA_TPCON

      SUBROUTINE DTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK,   &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTPCON

      SUBROUTINE ZTPCON( NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK,   &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTPCON

      END INTERFACE

      INTERFACE LA_TPTRS

      SUBROUTINE DTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DTPTRS

      SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZTPTRS

      MODULE PROCEDURE DTPTRS1
      MODULE PROCEDURE ZTPTRS1

      END INTERFACE

      INTERFACE LA_TRTRI

      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE DTRTRI

      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A( LDA, * )
      END SUBROUTINE ZTRTRI

      END INTERFACE

      INTERFACE LA_TRRFS

      SUBROUTINE DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, &
     &                   LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(IN) :: X(LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTRRFS

      SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X, &
     &                   LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTRRFS

      MODULE PROCEDURE DTRRFS1
      MODULE PROCEDURE ZTRRFS1

      END INTERFACE

      INTERFACE LA_TRCON

      SUBROUTINE DTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,      &
     &                   IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(IN) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DTRCON

      SUBROUTINE ZTRCON( NORM, UPLO, DIAG, N, A, LDA, RCOND, WORK,      &
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, NORM, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZTRCON

      END INTERFACE

      INTERFACE LA_TRTRS

      SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,    &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DTRTRS

      SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,    &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZTRTRS

      MODULE PROCEDURE DTRTRS1
      MODULE PROCEDURE ZTRTRS1

      END INTERFACE

      INTERFACE LA_SPTRI

      SUBROUTINE DSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSPTRI

      SUBROUTINE ZSPTRI( UPLO, N, AP, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSPTRI

      END INTERFACE

      INTERFACE LA_HPTRI

      SUBROUTINE ZHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHPTRI

      END INTERFACE

      INTERFACE LA_SPRFS

      SUBROUTINE DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,  &
     &                   FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: X(LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSPRFS

      SUBROUTINE ZSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,  &
     &                   FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSPRFS

      MODULE PROCEDURE DSPRFS1
      MODULE PROCEDURE ZSPRFS1

      END INTERFACE

      INTERFACE LA_HPRFS

      SUBROUTINE ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,  &
     &                   FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHPRFS

      MODULE PROCEDURE ZHPRFS1

      END INTERFACE

      INTERFACE LA_HPCON

      SUBROUTINE ZHPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV( * )
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHPCON

      END INTERFACE

      INTERFACE LA_SPCON

      SUBROUTINE DSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, IWORK,  &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV( * )
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSPCON

      SUBROUTINE ZSPCON( UPLO, N, AP, IPIV, ANORM, RCOND, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV( * )
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSPCON

      END INTERFACE

      INTERFACE LA_SPTRS

      SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DSPTRS

      SUBROUTINE ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZSPTRS

      MODULE PROCEDURE DSPTRS1
      MODULE PROCEDURE ZSPTRS1

      END INTERFACE

      INTERFACE LA_HPTRS

      SUBROUTINE ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZHPTRS

      MODULE PROCEDURE ZHPTRS1

      END INTERFACE

      INTERFACE LA_HPTRF

      SUBROUTINE ZHPTRF( UPLO, N, AP, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE ZHPTRF

      END INTERFACE

      INTERFACE LA_SPTRF

      SUBROUTINE DSPTRF( UPLO, N, AP, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE DSPTRF

      SUBROUTINE ZSPTRF( UPLO, N, AP, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE ZSPTRF

      END INTERFACE

      INTERFACE LA_SYTRI

      SUBROUTINE DSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: A( LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYTRI

      SUBROUTINE ZSYTRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSYTRI

      END INTERFACE

      INTERFACE LA_HETRI

      SUBROUTINE ZHETRI( UPLO, N, A, LDA, IPIV, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHETRI

      END INTERFACE

      INTERFACE LA_SYRFS

      SUBROUTINE DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, &
     &                   X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         INTEGER, INTENT(IN) :: IPIV(*) 
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B( LDB,*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYRFS

      SUBROUTINE ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, &
     &                   X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSYRFS

      MODULE PROCEDURE DSYRFS1
      MODULE PROCEDURE ZSYRFS1

      END INTERFACE

      INTERFACE LA_HERFS

      SUBROUTINE ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, &
     &                   X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHERFS

      MODULE PROCEDURE ZHERFS1

      END INTERFACE

      INTERFACE LA_SYCON

      SUBROUTINE DSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,     &
     &                   IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: A( LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYCON

      SUBROUTINE ZSYCON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,     &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSYCON

      END INTERFACE

      INTERFACE LA_HECON

      SUBROUTINE ZHECON( UPLO, N, A, LDA, IPIV, ANORM, RCOND, WORK,     &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHECON

      END INTERFACE

      INTERFACE LA_HETRS

      SUBROUTINE ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZHETRS

      MODULE PROCEDURE ZHETRS1

      END INTERFACE

      INTERFACE LA_SYTRS

      SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: A( LDA,*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DSYTRS

      SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZSYTRS

      MODULE PROCEDURE DSYTRS1
      MODULE PROCEDURE ZSYTRS1

      END INTERFACE

      INTERFACE LA_HETRF

      SUBROUTINE ZHETRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK( LWORK )
      END SUBROUTINE ZHETRF

      END INTERFACE

      INTERFACE LA_SYTRF

      SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         REAL(WP), INTENT(INOUT) :: A( LDA,*)
         REAL(WP), INTENT(OUT) :: WORK( LWORK )
      END SUBROUTINE DSYTRF

      SUBROUTINE ZSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A( LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK( LWORK )
      END SUBROUTINE ZSYTRF

      END INTERFACE

      INTERFACE LA_PTRFS

      SUBROUTINE DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR,   &
     &                   BERR, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*), DF(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: B( LDB,*), E(*), EF(*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPTRFS

      SUBROUTINE ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,   &
     &                   FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*), DF(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: B( LDB,*), E(*), EF(*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPTRFS

      MODULE PROCEDURE DPTRFS1
      MODULE PROCEDURE ZPTRFS1

      END INTERFACE

      INTERFACE LA_PTCON

      SUBROUTINE DPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM, D(*)
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         REAL(WP), INTENT(IN) :: E(*)
      END SUBROUTINE DPTCON

      SUBROUTINE ZPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM, D(*)
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: E(*)
      END SUBROUTINE ZPTCON

      END INTERFACE

      INTERFACE LA_PTTRS

      SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(IN) :: E(*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DPTTRS

      SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*)
         COMPLEX(WP), INTENT(IN) :: E(*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZPTTRS

      MODULE PROCEDURE DPTTRS1
      MODULE PROCEDURE ZPTTRS1

      END INTERFACE

      INTERFACE LA_PTTRF

      SUBROUTINE DPTTRF( N, D, E, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D( * )
         REAL(WP), INTENT(INOUT) :: E( * )
      END SUBROUTINE DPTTRF

      SUBROUTINE ZPTTRF( N, D, E, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D( * )
         COMPLEX(WP), INTENT(INOUT) :: E( * )
      END SUBROUTINE ZPTTRF

      END INTERFACE

      INTERFACE LA_PBEQU

      SUBROUTINE DPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
      END SUBROUTINE DPBEQU

      SUBROUTINE ZPBEQU( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
      END SUBROUTINE ZPBEQU

      END INTERFACE

      INTERFACE LA_PBRFS

      SUBROUTINE DPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B,    &
     &                   LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*), B( LDB,*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPBRFS

      SUBROUTINE ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B,    &
     &                   LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*),        &
     &                               B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPBRFS

      MODULE PROCEDURE DPBRFS1
      MODULE PROCEDURE ZPBRFS1

      END INTERFACE

      INTERFACE LA_PBCON

      SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,     &
     &                   IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPBCON

      SUBROUTINE ZPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,     &
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPBCON

      END INTERFACE

      INTERFACE LA_PBTRS

      SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DPBTRS

      SUBROUTINE ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZPBTRS

      MODULE PROCEDURE DPBTRS1
      MODULE PROCEDURE ZPBTRS1

      END INTERFACE

      INTERFACE LA_PBTRF

      SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB( LDAB,*)
      END SUBROUTINE DPBTRF

      SUBROUTINE ZPBTRF( UPLO, N, KD, AB, LDAB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB( LDAB,*)
      END SUBROUTINE ZPBTRF

      END INTERFACE

      INTERFACE LA_PPEQU

      SUBROUTINE DPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         REAL(WP), INTENT(IN) :: AP(*)
      END SUBROUTINE DPPEQU

      SUBROUTINE ZPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
      END SUBROUTINE ZPPEQU

      END INTERFACE

      INTERFACE LA_PPTRI

      SUBROUTINE DPPTRI( UPLO, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE DPPTRI

      SUBROUTINE ZPPTRI( UPLO, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE ZPPTRI

      END INTERFACE

      INTERFACE LA_PPRFS

      SUBROUTINE DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR,  &
     &                   BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: AFP(*), AP(*), B( LDB,*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPPRFS

      SUBROUTINE ZPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR,  &
     &                   BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPPRFS

      MODULE PROCEDURE DPPRFS1
      MODULE PROCEDURE ZPPRFS1

      END INTERFACE

      INTERFACE LA_PPCON

      SUBROUTINE DPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPPCON

      SUBROUTINE ZPPCON( UPLO, N, AP, ANORM, RCOND, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPPCON

      END INTERFACE

      INTERFACE LA_PPTRS

      SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DPPTRS

      SUBROUTINE ZPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZPPTRS

      MODULE PROCEDURE DPPTRS1
      MODULE PROCEDURE ZPPTRS1

      END INTERFACE

      INTERFACE LA_PPTRF

      SUBROUTINE DPPTRF( UPLO, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE DPPTRF

      SUBROUTINE ZPPTRF( UPLO, N, AP, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
      END SUBROUTINE ZPPTRF

      END INTERFACE

      INTERFACE LA_POEQU

      SUBROUTINE DPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         REAL(WP), INTENT(IN) :: A( LDA,*)
      END SUBROUTINE DPOEQU

      SUBROUTINE ZPOEQU( N, A, LDA, S, SCOND, AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, SCOND, S(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
      END SUBROUTINE ZPOEQU

      END INTERFACE

      INTERFACE LA_POTRI

      SUBROUTINE DPOTRI( UPLO, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A( LDA,*)
      END SUBROUTINE DPOTRI

      SUBROUTINE ZPOTRI( UPLO, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A( LDA,*)
      END SUBROUTINE ZPOTRI

      END INTERFACE

      INTERFACE LA_PORFS

      SUBROUTINE DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,    &
     &                   LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: A( LDA,*), AF( LDAF,*), B( LDB,*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DPORFS

      SUBROUTINE ZPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,    &
     &                   LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*), AF( LDAF,*), B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZPORFS

      MODULE PROCEDURE DPORFS1
      MODULE PROCEDURE ZPORFS1

      END INTERFACE

      INTERFACE LA_POTRS

      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A( LDA,*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DPOTRS

      SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZPOTRS

      MODULE PROCEDURE DPOTRS1
      MODULE PROCEDURE ZPOTRS1

      END INTERFACE

      INTERFACE LA_GTRFS

      SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,  &
     &                   IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: B( LDB,*), D(*), DF(*), DL(*), DLF(*), &
     &                           DU(*), DU2(*), DUF(*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGTRFS

      SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,  &
     &                   IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: B( LDB,*), D(*), DF(*), DL(*),      &
     &                              DLF(*), DU(*), DU2(*), DUF(*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZGTRFS

      MODULE PROCEDURE DGTRFS1
      MODULE PROCEDURE ZGTRFS1

      END INTERFACE

      INTERFACE LA_GTCON

      SUBROUTINE DGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,   &
     &                   WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGTCON

      SUBROUTINE ZGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,   &
     &                   WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZGTCON

      END INTERFACE

      INTERFACE LA_GTTRS

      SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,  &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         REAL(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE DGTTRS

      SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,  &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
      END SUBROUTINE ZGTTRS

      MODULE PROCEDURE DGTTRS1
      MODULE PROCEDURE ZGTTRS1

      END INTERFACE

      INTERFACE LA_GTTRF

      SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: D(*), DL(*), DU(*)
         REAL(WP), INTENT(OUT) :: DU2(*)
      END SUBROUTINE DGTTRF

      SUBROUTINE ZGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: D(*), DL(*), DU(*)
         COMPLEX(WP), INTENT(OUT) :: DU2(*)
      END SUBROUTINE ZGTTRF

      END INTERFACE

      INTERFACE LA_GBEQU

      SUBROUTINE DGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,  &
     &                   AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, COLCND, ROWCND
         REAL(WP), INTENT(OUT) :: C(*), R(*)
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
      END SUBROUTINE DGBEQU

      SUBROUTINE ZGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,  &
     &                   AMAX, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, COLCND, ROWCND
         REAL(WP), INTENT(OUT) :: C(*), R(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
      END SUBROUTINE ZGBEQU

      END INTERFACE

      INTERFACE LA_GBRFS

      SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,  &
     &                   IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
         REAL(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*), B( LDB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(INOUT) :: X( LDX,*)
      END SUBROUTINE DGBRFS

      SUBROUTINE ZGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,  &
     &                   IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*),         &
     &                              B( LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
      END SUBROUTINE ZGBRFS

      MODULE PROCEDURE DGBRFS1
      MODULE PROCEDURE ZGBRFS1
      END INTERFACE

      INTERFACE LA_GBCON

      SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, &
     &                   WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: KL, KU, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV( * )
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(WP), INTENT(IN) :: AB( LDAB, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DGBCON

      SUBROUTINE ZGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND, &
     &                   WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: KL, KU, LDAB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(IN) :: IPIV( * )
         REAL(WP), INTENT(OUT) :: RWORK( * )
         COMPLEX(WP), INTENT(IN) :: AB( LDAB, * )
         COMPLEX(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZGBCON

      END INTERFACE

      INTERFACE LA_GBTRS

      SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DGBTRS

      SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZGBTRS

      MODULE PROCEDURE DGBTRS1
      MODULE PROCEDURE ZGBTRS1
      END INTERFACE

      INTERFACE LA_GBTRF

      SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*)
      END SUBROUTINE DGBTRF

      SUBROUTINE ZGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*)
      END SUBROUTINE ZGBTRF

      END INTERFACE

      INTERFACE

      FUNCTION DLAMCH( CMACH )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLAMCH
         CHARACTER(LEN=1), INTENT(IN) :: CMACH
      END FUNCTION DLAMCH

      END INTERFACE

      INTERFACE LA_GGSVD

       SUBROUTINE DGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,   &
     &                    LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,     &
     &                    WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBV, JOBQ
         INTEGER, INTENT(IN) :: M, N, P, LDA, LDB, LDU, LDV, LDQ
         INTEGER, INTENT(OUT) :: INFO, K, L, IWORK(*)
         REAL(WP), INTENT(OUT) :: ALPHA(*), BETA(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: U(LDU,*), V(LDV,*), Q(LDQ,*), WORK(*)
      END SUBROUTINE DGGSVD

       SUBROUTINE ZGGSVD( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,   &
     &                    LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,     &
     &                    WORK, RWORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBV, JOBQ
         INTEGER, INTENT(IN) :: M, N, P, LDA, LDB, LDU, LDV, LDQ
         INTEGER, INTENT(OUT) :: INFO, K, L, IWORK(*)
         REAL(WP), INTENT(OUT) :: ALPHA(*), BETA(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: U(LDU,*), V(LDV,*), Q(LDQ,*),      &
     &                               WORK(*)
      END SUBROUTINE ZGGSVD

       END INTERFACE

      INTERFACE LA_GEGV

       SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,       &
     &                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VL(LDVL,*), VR(LDVR,*), WORK(*)
      END SUBROUTINE DGEGV

       SUBROUTINE ZGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,  &
     &                   VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), VL(LDVL,*),     &
     &                               VR(LDVR,*), WORK(*)
      END SUBROUTINE ZGEGV

       END INTERFACE

      INTERFACE LA_GEGS

       SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR,     &
     &                   ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK,    &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VSL(LDVSL,*), VSR(LDVSR,*), WORK(*)
      END SUBROUTINE DGEGS

       SUBROUTINE ZGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA,&
     &                   VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK,    &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), VSL(LDVSL,*),   &
     &                               VSR(LDVSR,*), WORK(*)
      END SUBROUTINE ZGEGS

       END INTERFACE

        INTERFACE LA_SBGVX
	

       SUBROUTINE DSBGVX( JOBZ, RANGE, UPLO, N, KAB, KBB, AB, LDAB, BB, &
     &                    LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,&
     &                    LDZ, WORK, IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: N, IL, IU, KAB, KBB, LDAB, LDBB, LDQ, LDZ
         INTEGER, INTENT(OUT) :: M
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
         REAL(WP), INTENT(OUT) :: WORK(*), Q(LDQ,*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
         INTEGER, INTENT(IN) :: IFAIL(*)
        END SUBROUTINE DSBGVX
	
        END INTERFACE
	    
       INTERFACE LA_HBGVX
		    
       SUBROUTINE ZHBGVX( JOBZ, RANGE, UPLO, N, KAB, KBB, AB, LDAB, BB, &
     &                    LDBB, Q, LDQ, VL, VU, IL, IU, ABSTOL, M, W, Z,&
     &                    LDZ, WORK, RWORK, IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: N, IL, IU, KAB, KBB, LDAB, LDBB, LDQ, LDZ
         INTEGER, INTENT(OUT) :: M
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), Q(LDQ,*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: RWORK (*)
         INTEGER, INTENT(IN) :: IFAIL(*)
        END SUBROUTINE ZHBGVX
	
        END INTERFACE

        INTERFACE LA_SBGVD
	

        SUBROUTINE DSBGVD( JOBZ, UPLO, N, KAB, KBB, AB, LDAB, BB, LDBB, &
     &                     W, Z, LDZ, WORK, LWORK, IWORK, LIWORK,       &
     &                     INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO 
           INTEGER, INTENT(IN) :: N, KAB, KBB, LDAB, LDBB, LDZ
           INTEGER, INTENT(IN) :: LWORK, LIWORK
           INTEGER, INTENT(OUT) :: INFO, IWORK(*)
           REAL(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
           REAL(WP), INTENT(OUT) :: W(*)
           REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
	 END SUBROUTINE DSBGVD
	 
         END INTERFACE
		
         INTERFACE LA_HBGVD
	 
       SUBROUTINE ZHBGVD( JOBZ, UPLO, N, KAB, KBB, AB, LDAB, BB, LDBB,  &
     &                    W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, &
     &                    LIWORK, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
           INTEGER, INTENT(IN) :: N, KAB, KBB, LDAB, LDBB, LDZ
           INTEGER, INTENT(IN) :: LWORK, LRWORK, LIWORK
           INTEGER, INTENT(OUT) :: INFO, IWORK(*)
           COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
           REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
           COMPLEX(WP),INTENT(OUT) :: Z(LDZ,*), WORK(*)
          END SUBROUTINE ZHBGVD
	  
        
         END INTERFACE

      INTERFACE LA_SBGV

       SUBROUTINE DSBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,  &
     &                   Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: KA, KB, LDAB, LDBB, N, LDZ
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
         REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE DSBGV

       END INTERFACE

      INTERFACE LA_HBGV

       SUBROUTINE ZHBGV( JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W,  &
     &                   Z, LDZ, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: KA, KB, LDAB, LDBB, N, LDZ
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), BB(LDBB,*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHBGV

       END INTERFACE

       INTERFACE LA_SPGVX


       SUBROUTINE DSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,  &
     &                    IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK,    &
     &                    IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, IL, IU, LDZ
         INTEGER, INTENT(OUT) :: M
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*), BP(*)
         REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
         INTEGER, INTENT(IN) :: IFAIL(*)
        END SUBROUTINE DSPGVX

        END INTERFACE
		
        INTERFACE LA_HPGVX
			
       SUBROUTINE ZHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU,  &
     &                    IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK,    &
     &                    IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, IL, IU, LDZ
         INTEGER, INTENT(OUT) :: M
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*), BP(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         INTEGER, INTENT(IN) :: IFAIL(*)
        END SUBROUTINE ZHPGVX
	

        END INTERFACE

       INTERFACE LA_SPGVD
       

        SUBROUTINE DSPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ,     &
     &                     WORK, LWORK, IWORK, LIWORK, INFO )
        USE LA_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
        INTEGER, INTENT(IN) :: ITYPE, N, LDZ
        INTEGER, INTENT(IN) :: LWORK, LIWORK
        INTEGER, INTENT(OUT) :: INFO, IWORK(*)
        REAL(WP), INTENT(INOUT) :: AP(*), BP(*)
        REAL(WP), INTENT(OUT) :: W(*)
        REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
       END SUBROUTINE DSPGVD

        END INTERFACE
		
        INTERFACE LA_HPGVD 

       SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK,&
     &                    LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
        USE LA_PRECISION, ONLY: WP => DP
        CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
        INTEGER, INTENT(IN) :: ITYPE, N, LDZ
        INTEGER, INTENT(IN) :: LWORK, LRWORK, LIWORK
        INTEGER, INTENT(OUT) :: INFO, IWORK(*)
        COMPLEX(WP), INTENT(INOUT) :: AP(*), BP(*)
        REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
        COMPLEX(WP), INTENT(OUT):: Z(LDZ,*), WORK(*)
       END SUBROUTINE ZHPGVD

        END INTERFACE

      INTERFACE LA_SPGV

       SUBROUTINE DSPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, LDZ
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(INOUT) :: AP(*), BP(*)
         REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPGV

       END INTERFACE

      INTERFACE LA_HPGV

       SUBROUTINE ZHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, &
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, LDZ
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), BP(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHPGV

       END INTERFACE

      INTERFACE LA_GESVD

       SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT,     &
     &                    LDVT, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBVT
         INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: S(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
      END SUBROUTINE DGESVD

       SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT,     &
     &                    LDVT, WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBU, JOBVT
         INTEGER, INTENT(IN) :: M, N, LDA, LDU, LDVT, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: S(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: U(LDU,*), VT(LDVT,*), WORK(*)
      END SUBROUTINE ZGESVD

       END INTERFACE

      INTERFACE LA_GEEVX

       SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR,   &
     &                    WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE,      &
     &                    ABNRM, RCONDE, RCONDV, WORK, LWORK, IWORK,    &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: BALANC, JOBVL, JOBVR, SENSE
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO, ILO, IHI, IWORK(*)
         REAL(WP), INTENT(OUT) :: ABNRM
         REAL(WP), INTENT(OUT) :: SCALE(*), RCONDE(*), RCONDV(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), WR(*), WI(*), &
     &                            WORK(*)
      END SUBROUTINE DGEEVX

       SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL,&
     &                    LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,       &
     &                    RCONDE, RCONDV, WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: BALANC, JOBVL, JOBVR, SENSE
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO, ILO, IHI
         REAL(WP), INTENT(OUT) :: ABNRM
         REAL(WP), INTENT(OUT) :: SCALE(*), RCONDE(*), RCONDV(*),       &
     &                            RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), W(*),      &
     &                               WORK(*)
      END SUBROUTINE ZGEEVX

       END INTERFACE

        INTERFACE LA_GGEVX

       SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B,    &
     &                    LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR,&
     &                    ILO, IHI, LSCALE, RSCALE, ABNRM, BBNRM,       &
     &                    RCONDE, RCONDV, WORK, LWORK, IWORK, BWORK,    &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: BALANC, JOBVL, JOBVR, SENSE
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT):: ILO, IHI
         REAL(WP), INTENT(OUT) :: ABNRM, BBNRM
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VL(LDVL,*), VR(LDVR,*), WORK(*),      &
     &                            LSCALE(*), RSCALE(*), RCONDE(*),      &
     &                            RCONDV(*)
         INTEGER :: IWORK(*)
         LOGICAL :: BWORK(*) 
        END SUBROUTINE DGGEVX

       SUBROUTINE ZGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B,    &
     &                    LDB, ALPHA, BETA, VL, LDVL, VR, LDVR, ILO,    &
     &                    IHI, LSCALE, RSCALE, ABNRM, BBNRM, RCONDE,    &
     &                    RCONDV, WORK, LWORK, RWORK, IWORK, BWORK,     &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: BALANC, JOBVL, JOBVR, SENSE
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT):: ILO, IHI
         REAL(WP), INTENT(OUT) :: ABNRM, BBNRM
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), VL(LDVL,*),     &
     &                               VR(LDVR,*), WORK(*)
         REAL(WP), INTENT(OUT) :: LSCALE(*), RSCALE(*), RCONDE(*), RCONDV(*)
         INTEGER :: IWORK(*)
         LOGICAL :: BWORK(*)
         REAL(WP) :: RWORK(*)
        END SUBROUTINE ZGGEVX

        END INTERFACE

        INTERFACE LA_GGEV
      
       SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,       &
     &                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VL(LDVL,*), VR(LDVR,*), WORK(*)
       END SUBROUTINE DGGEV

       SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,  &
     &                   VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
       INTEGER, INTENT(IN) :: LDA, LDB, N, LDVL, LDVR, LWORK
       INTEGER, INTENT(OUT) :: INFO
       COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) ::  ALPHA(*), BETA(*), VL(LDVL,*),    &
     &                                VR(LDVR,*), WORK(*)
       REAL(WP) :: RWORK(*)
      END SUBROUTINE ZGGEV
      
       END INTERFACE

      INTERFACE LA_GEEV

       SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
     &                   LDVR, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), WR(*), WI(*), &
     &                            WORK(*)
      END SUBROUTINE DGEEV

       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,&
     &                   WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVL, JOBVR
         INTEGER, INTENT(IN) :: N, LDA, LDVL, LDVR, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: VL(LDVL,*), VR(LDVR,*), W(*),      &
     &                               WORK(*)
      END SUBROUTINE ZGEEV

       END INTERFACE

      INTERFACE LA_GEESX

       SUBROUTINE DGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM,  &
     &                    WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK,&
     &                    IWORK, LIWORK, BWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTERFACE
            LOGICAL FUNCTION SELECT(WR, WI)
               USE LA_PRECISION, ONLY: WP => DP
               REAL(WP), INTENT(IN) :: WR, WI
            END FUNCTION SELECT
         END INTERFACE
         CHARACTER(LEN=1), INTENT(IN) :: JOBVS, SORT, SENSE
         INTEGER, INTENT(IN) :: N, LDA, LDVS, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, SDIM, IWORK(*)
         LOGICAL, INTENT(OUT) :: BWORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: RCONDV, RCONDE
         REAL(WP), INTENT(OUT) :: VS(LDVS,*), WR(*), WI(*), WORK(*)
         OPTIONAL :: SELECT
      END SUBROUTINE DGEESX

       SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM,  &
     &                    W, VS, LDVS, RCONDE, RCONDV, WORK, LWORK,     &
     &                    RWORK, BWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTERFACE
            LOGICAL FUNCTION SELECT( W )
               USE LA_PRECISION, ONLY: WP => DP
               COMPLEX(WP), INTENT(IN) :: W
            END FUNCTION SELECT
         END INTERFACE
         CHARACTER(LEN=1), INTENT(IN) :: JOBVS, SORT, SENSE
         INTEGER, INTENT(IN) :: N, LDA, LDVS, LWORK
         INTEGER, INTENT(OUT) :: INFO, SDIM
         LOGICAL, INTENT(OUT) :: BWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: RCONDV, RCONDE
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: VS(LDVS,*), W(*), WORK(*)
         OPTIONAL :: SELECT
      END SUBROUTINE ZGEESX

       END INTERFACE
       
       INTERFACE LA_GGESX

       SUBROUTINE DGGESX( JOBVSL, JOBVSR, SORT, DELCTG, SENSE, N, A,    &
     &                    LDA, B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, &
     &                    LDVSL, VSR, LDVSR, RCONDE, RCONDV, WORK,      &
     &                    LWORK, IWORK, LIWORK, BWORK, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR, SORT, SENSE
          INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK, LIWORK
          INTEGER, INTENT(INOUT) :: INFO
          INTEGER, INTENT(OUT) :: SDIM
          REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VSL(LDVSL,*), VSR(LDVSR,*), WORK(*)
          REAL(WP), INTENT(OUT) :: RCONDE(2), RCONDV(2)
          LOGICAL :: BWORK(*)
          INTEGER :: IWORK (*)
          INTERFACE
           LOGICAL FUNCTION DELCTG(ALPHAR, ALPHAI, BETA)
              USE LA_PRECISION, ONLY: WP => DP
              REAL(WP), INTENT(IN) :: ALPHAR, ALPHAI, BETA
           END FUNCTION DELCTG
          END INTERFACE
          OPTIONAL :: DELCTG
         END SUBROUTINE DGGESX
	 
       SUBROUTINE ZGGESX( JOBVSL, JOBVSR, SORT, DELCTG, SENSE, N, A,    &
     &                    LDA, B, LDB, SDIM, ALPHA, BETA, VSL, LDVSL,   &
     &                    VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK,      &
     &                    RWORK, IWORK, LIWORK, BWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR, SORT, SENSE
         INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK, LIWORK
         INTEGER, INTENT(INOUT) :: INFO
         INTEGER, INTENT(OUT) :: SDIM
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), VSL(LDVSL,*),   &
     &                               VSR(LDVSR,*), WORK(*)
         REAL(WP), INTENT(OUT) :: RCONDE(2), RCONDV(2)
         LOGICAL :: BWORK(*)
         INTEGER :: IWORK (*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         INTERFACE
          LOGICAL FUNCTION DELCTG(ALPHA, BETA)
            USE LA_PRECISION, ONLY: WP => DP
            COMPLEX(WP), INTENT(IN) :: ALPHA, BETA
          END FUNCTION DELCTG
        END INTERFACE
        OPTIONAL :: DELCTG
      END SUBROUTINE ZGGESX

      END INTERFACE 

       
      INTERFACE LA_GGES
      
       SUBROUTINE DGGES( JOBVSL, JOBVSR, SORT, DELCTG, N, A, LDA, B,    &
     &                   LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL,   &
     &                   VSR, LDVSR, WORK, LWORK, BWORK, INFO )
       USE LA_PRECISION, ONLY: WP => DP
       CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR, SORT
       INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK
       INTEGER, INTENT(INOUT) :: INFO
       INTEGER, INTENT(OUT) :: SDIM
       REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: ALPHAR(*), ALPHAI(*), BETA(*),        &
     &                            VSL(LDVSL,*), VSR(LDVSR,*), WORK(*)
       LOGICAL :: BWORK(*)
       INTERFACE
         LOGICAL FUNCTION DELCTG(ALPHAR, ALPHAI, BETA)
          USE LA_PRECISION, ONLY: WP => DP
          REAL(WP), INTENT(IN) :: ALPHAR, ALPHAI, BETA
         END FUNCTION DELCTG
       END INTERFACE
       OPTIONAL :: DELCTG
      END SUBROUTINE DGGES
      
       SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, DELCTG, N, A, LDA, B,    &
     &                   LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR,&
     &                   WORK, LWORK, RWORK, BWORK, INFO )
       USE LA_PRECISION, ONLY : WP => DP
       CHARACTER(LEN=1), INTENT(IN) :: JOBVSL, JOBVSR, SORT
       INTEGER, INTENT(IN) :: LDA, LDB, N, LDVSL, LDVSR, LWORK
       INTEGER, INTENT(INOUT) :: INFO
       INTEGER, INTENT(OUT) :: SDIM
       COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: ALPHA(*), BETA(*), VSL(LDVSL,*),   &
     &                               VSR(LDVSR,*), WORK(*)
       LOGICAL :: BWORK(*)
       REAL(WP) :: RWORK(*)
       INTERFACE
         LOGICAL FUNCTION DELCTG( ALPHA, BETA)
           USE LA_PRECISION, ONLY: WP => DP
           COMPLEX(WP), INTENT(IN) :: ALPHA, BETA
         END FUNCTION DELCTG
       END INTERFACE
       OPTIONAL :: DELCTG
      END SUBROUTINE ZGGES
      
        END INTERFACE

      INTERFACE LA_GEES

       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,  &
     &                   VS, LDVS, WORK, LWORK, BWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTERFACE
            LOGICAL FUNCTION SELECT(WR, WI)
               USE LA_PRECISION, ONLY: WP => DP
               REAL(WP), INTENT(IN) :: WR, WI
            END FUNCTION SELECT
         END INTERFACE
         CHARACTER(LEN=1), INTENT(IN) :: JOBVS, SORT
         INTEGER, INTENT(IN) :: N, LDA, LDVS, LWORK
         INTEGER, INTENT(OUT) :: INFO, SDIM
         LOGICAL, INTENT(OUT) :: BWORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: VS(LDVS,*), WR(*), WI(*), WORK(*)
         OPTIONAL :: SELECT
      END SUBROUTINE DGEES

       SUBROUTINE ZGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, W, VS,   &
     &                   LDVS, WORK, LWORK, RWORK, BWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTERFACE
            LOGICAL FUNCTION SELECT( W )
               USE LA_PRECISION, ONLY: WP => DP
               COMPLEX(WP), INTENT(IN) :: W
            END FUNCTION SELECT
         END INTERFACE
         CHARACTER(LEN=1), INTENT(IN) :: JOBVS, SORT
         INTEGER, INTENT(IN) :: N, LDA, LDVS, LWORK
         INTEGER, INTENT(OUT) :: INFO, SDIM
         LOGICAL, INTENT(OUT) :: BWORK(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: VS(LDVS,*), W(*), WORK(*)
         OPTIONAL :: SELECT
      END SUBROUTINE ZGEES

       END INTERFACE

        INTERFACE LA_STEVR

       SUBROUTINE DSTEVR( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, &
     &                    M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,     &
     &                    LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE
         INTEGER, INTENT(IN) :: N, IL, IU, LDZ,  LWORK, LIWORK
         INTEGER, INTENT(OUT) :: M
         INTEGER, INTENT(OUT), TARGET :: ISUPPZ(*)
         REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: WORK(*), W(*)
         REAL(WP), INTENT(OUT), TARGET :: Z(LDZ,*)
       END SUBROUTINE DSTEVR

       END INTERFACE

      INTERFACE LA_STEVX

       SUBROUTINE DSTEVX( JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, &
     &                    M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSTEVX

       END INTERFACE

      INTERFACE LA_STEVD

       SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK,    &
     &                    LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ
         INTEGER, INTENT(IN) :: LDZ, N, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE DSTEVD

       END INTERFACE

      INTERFACE LA_STEV

       SUBROUTINE DSTEV( JOBZ, N, D, E, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE DSTEV

       END INTERFACE

      INTERFACE LA_SBEVX

       SUBROUTINE DSBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ,   &
     &                    VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,   &
     &                    IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU, LDQ, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: Q(LDQ,*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSBEVX

       END INTERFACE

      INTERFACE LA_HBEVX

       SUBROUTINE ZHBEVX( JOBZ, RANGE, UPLO, N, KD, AB, LDAB, Q, LDQ,   &
     &                    VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,   &
     &                    RWORK, IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU, LDQ, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Q(LDQ,*), Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHBEVX

       END INTERFACE

      INTERFACE LA_SBEVD

       SUBROUTINE DSBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, &
     &                    LWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: N, KD, LDAB, LDZ, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE DSBEVD

       END INTERFACE

      INTERFACE LA_HBEVD

       SUBROUTINE ZHBEVD( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK, &
     &                    LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: N, KD, LDAB, LDZ, LWORK, LIWORK, LRWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHBEVD

       END INTERFACE

      INTERFACE LA_SBEV

       SUBROUTINE DSBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,  &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: N, KD, LDAB, LDZ
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSBEV

       END INTERFACE

      INTERFACE LA_HBEV

       SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,  &
     &                   RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: N, KD, LDAB, LDZ
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHBEV

       END INTERFACE

      INTERFACE LA_SPEVX

       SUBROUTINE DSPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,     &
     &                    ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL,     &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEVX

       END INTERFACE

      INTERFACE LA_HPEVX

       SUBROUTINE ZHPEVX( JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU,     &
     &                    ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK,     &
     &                    IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO, RANGE
         INTEGER, INTENT(IN) :: LDZ, N, IL, IU
         INTEGER, INTENT(OUT) :: INFO, IWORK(*), M, IFAIL(*)
         REAL(WP), INTENT(IN) :: VL, VU, ABSTOL
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHPEVX

       END INTERFACE

      INTERFACE LA_SPEVD

       SUBROUTINE DSPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,    &
     &                    IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEVD

       END INTERFACE

      INTERFACE LA_HPEVD

       SUBROUTINE ZHPEVD( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK,    &
     &                    RWORK, LRWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N, LWORK, LIWORK, LRWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHPEVD

       END INTERFACE

      INTERFACE LA_SPEV

       SUBROUTINE DSPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), Z(LDZ,*), WORK(*)
      END SUBROUTINE DSPEV

       END INTERFACE

      INTERFACE LA_HPEV

       SUBROUTINE ZHPEV( JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK,     &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDZ, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: Z(LDZ,*), WORK(*)
      END SUBROUTINE ZHPEV

       END INTERFACE

      INTERFACE LA_GGGLM

       SUBROUTINE DGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: P, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), D(*)
         REAL(WP), INTENT(OUT) :: WORK(*), X(*), Y(*)
      END SUBROUTINE DGGGLM

       SUBROUTINE ZGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: P, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), D(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), X(*), Y(*)
      END SUBROUTINE ZGGGLM

       END INTERFACE

      INTERFACE LA_GGLSE

       SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: P, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), C(*), D(*)
         REAL(WP), INTENT(OUT) :: WORK(*), X(*)
      END SUBROUTINE DGGLSE

       SUBROUTINE ZGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: P, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*), C(*), D(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), X(*)
      END SUBROUTINE ZGGLSE

       END INTERFACE

       INTERFACE LA_GELSY

       SUBROUTINE DGELSY(  M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,     &
     &                     RANK, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) ::  WORK(*)
       END SUBROUTINE DGELSY

       SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,&
     &                    WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) ::  WORK(*)
         REAL(WP) :: RWORK(*)
       END SUBROUTINE ZGELSY

        MODULE PROCEDURE DGELSY1
        MODULE PROCEDURE ZGELSY1
			
        END INTERFACE 

        INTERFACE LA_GELSD

       SUBROUTINE DGELSD(  M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, IWORK, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
          INTEGER, INTENT(OUT) :: INFO, RANK
          REAL(WP), INTENT(IN) :: RCOND
          REAL(WP), INTENT(INOUT) ::  A(LDA,*), B(LDB,*)
          REAL(WP), INTENT(OUT) :: S(*)
          REAL(WP), INTENT(OUT) :: WORK(*)
          INTEGER :: IWORK(*) 
        END SUBROUTINE DGELSD

       SUBROUTINE ZGELSD(  M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, RWORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: S(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         REAL(WP) :: RWORK(*)
         INTEGER :: IWORK(*)
       END SUBROUTINE ZGELSD

          MODULE PROCEDURE DGELSD1
          MODULE PROCEDURE ZGELSD1
		  
       END INTERFACE

      INTERFACE LA_GELSX

       SUBROUTINE DGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,&
     &                    WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGELSX

       SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,&
     &                    WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZGELSX

      MODULE PROCEDURE DGELSX1
      MODULE PROCEDURE ZGELSX1

       END INTERFACE

      INTERFACE LA_GELSS

       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,   &
     &                    WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: S(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGELSS

       SUBROUTINE ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,   &
     &                    WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: S(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZGELSS

      MODULE PROCEDURE DGELSS1
      MODULE PROCEDURE ZGELSS1

       END INTERFACE

      INTERFACE LA_GELS

       SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DGELS

       SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,&
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZGELS

      MODULE PROCEDURE DGELS1
      MODULE PROCEDURE ZGELS1

       END INTERFACE

      INTERFACE LA_SPSV

       SUBROUTINE DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
      END SUBROUTINE DSPSV

       SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
      END SUBROUTINE ZSPSV

      MODULE PROCEDURE DSPSV1
      MODULE PROCEDURE ZSPSV1

       END INTERFACE

      INTERFACE LA_HPSV

       SUBROUTINE ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
      END SUBROUTINE ZHPSV

      MODULE PROCEDURE ZHPSV1

       END INTERFACE

      INTERFACE LA_SYSV

       SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYSV

       SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZSYSV

      MODULE PROCEDURE DSYSV1
      MODULE PROCEDURE ZSYSV1

       END INTERFACE

      INTERFACE LA_HESV

       SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,     &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHESV

      MODULE PROCEDURE ZHESV1

       END INTERFACE

      INTERFACE LA_PTSV

       SUBROUTINE DPTSV( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*)
         REAL(WP), INTENT(INOUT) :: E(*), B(LDB,*)
      END SUBROUTINE DPTSV

       SUBROUTINE ZPTSV( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*)
         COMPLEX(WP), INTENT(INOUT) :: E(*), B(LDB,*)
      END SUBROUTINE ZPTSV

      MODULE PROCEDURE DPTSV1
      MODULE PROCEDURE ZPTSV1

       END INTERFACE

      INTERFACE LA_PBSV

       SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
      END SUBROUTINE DPBSV

       SUBROUTINE ZPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
      END SUBROUTINE ZPBSV

      MODULE PROCEDURE DPBSV1
      MODULE PROCEDURE ZPBSV1

       END INTERFACE

      INTERFACE LA_PPSV

       SUBROUTINE DPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
      END SUBROUTINE DPPSV

       SUBROUTINE ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
      END SUBROUTINE ZPPSV

      MODULE PROCEDURE DPPSV1
      MODULE PROCEDURE ZPPSV1

       END INTERFACE

      INTERFACE LA_POSV

       SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DPOSV

       SUBROUTINE ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE ZPOSV

      MODULE PROCEDURE DPOSV1
      MODULE PROCEDURE ZPOSV1

       END INTERFACE

      INTERFACE LA_GTSV

       SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(LDB,*)
      END SUBROUTINE DGTSV

       SUBROUTINE ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(LDB,*)
      END SUBROUTINE ZGTSV

      MODULE PROCEDURE DGTSV1
      MODULE PROCEDURE ZGTSV1

       END INTERFACE

      INTERFACE LA_GBSV

       SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
      END SUBROUTINE DGBSV

       SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
      END SUBROUTINE ZGBSV

      MODULE PROCEDURE DGBSV1
      MODULE PROCEDURE ZGBSV1

       END INTERFACE

      INTERFACE LA_GESV

       SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE DGESV

       SUBROUTINE ZGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
      END SUBROUTINE ZGESV

      MODULE PROCEDURE DGESV1
      MODULE PROCEDURE ZGESV1

       END INTERFACE

      INTERFACE LA_SPSVX

       SUBROUTINE DSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X,  &
     &                    LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(IN) :: A(*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: AF(*)
      END SUBROUTINE DSPSVX

       SUBROUTINE ZSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X,  &
     &                    LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: AF(*)
      END SUBROUTINE ZSPSVX

      MODULE PROCEDURE DSPSVX1
      MODULE PROCEDURE ZSPSVX1

       END INTERFACE

      INTERFACE LA_HPSVX

       SUBROUTINE ZHPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X,  &
     &                    LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: AF(*)
      END SUBROUTINE ZHPSVX

      MODULE PROCEDURE ZHPSVX1

       END INTERFACE

      INTERFACE LA_SYSVX

       SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,  &
     &                    B, LDB, X, LDX, RCOND, FERR, BERR, WORK,      &
     &                    LWORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: AF(LDAF,*)
      END SUBROUTINE DSYSVX

       SUBROUTINE ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,  &
     &                    B, LDB, X, LDX, RCOND, FERR, BERR, WORK,      &
     &                    LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
      END SUBROUTINE ZSYSVX

      MODULE PROCEDURE DSYSVX1
      MODULE PROCEDURE ZSYSVX1

       END INTERFACE

      INTERFACE LA_HESVX

       SUBROUTINE ZHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,  &
     &                    B, LDB, X, LDX, RCOND, FERR, BERR, WORK,      &
     &                    LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
      END SUBROUTINE ZHESVX

      MODULE PROCEDURE ZHESVX1

       END INTERFACE

      INTERFACE LA_PTSVX

       SUBROUTINE DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,  &
     &                    RCOND, FERR, BERR, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(IN) :: E(*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: DF(*)
         REAL(WP), INTENT(INOUT) :: EF(*)
      END SUBROUTINE DPTSVX

       SUBROUTINE ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,  &
     &                    RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(IN) :: D(*)
         COMPLEX(WP), INTENT(IN) :: E(*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: DF(*)
         COMPLEX(WP), INTENT(INOUT) :: EF(*)
      END SUBROUTINE ZPTSVX

      MODULE PROCEDURE DPTSVX1
      MODULE PROCEDURE ZPTSVX1

       END INTERFACE

      INTERFACE LA_PBSVX

       SUBROUTINE DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB,&
     &                    EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR,  &
     &                    WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*), FERR(*), BERR(*),  &
     &                            RCOND
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*), B(LDB,*), &
     &                              S(*)
      END SUBROUTINE DPBSVX

       SUBROUTINE ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB,&
     &                    EQUED, S, B, LDB, X, LDX, RCOND, FERR, BERR,  &
     &                    WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RCOND, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*),        &
     &                                 B(LDB,*)
      END SUBROUTINE ZPBSVX

      MODULE PROCEDURE DPBSVX1
      MODULE PROCEDURE ZPBSVX1

       END INTERFACE

      INTERFACE LA_PPSVX

       SUBROUTINE DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,    &
     &                    LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK,  &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: AP(*), AFP(*), B(LDB,*), S(*)
      END SUBROUTINE DPPSVX

       SUBROUTINE ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,    &
     &                    LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK,  &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), AFP(*), B(LDB,*)
      END SUBROUTINE ZPPSVX

      MODULE PROCEDURE DPPSVX1
      MODULE PROCEDURE ZPPSVX1

       END INTERFACE

      INTERFACE LA_POSVX

       SUBROUTINE DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, &
     &                    S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,   &
     &                    IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE DPOSVX

       SUBROUTINE ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, &
     &                    S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,   &
     &                    RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE ZPOSVX

      MODULE PROCEDURE DPOSVX1
      MODULE PROCEDURE ZPOSVX1

       END INTERFACE

      INTERFACE LA_GTSVX

       SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,&
     &                    DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, &
     &                    WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), X(LDX,*), WORK(*)
         REAL(WP), INTENT(IN) :: B(LDB,*), DL(*), D(*), DU(*) 
         REAL(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
      END SUBROUTINE DGTSVX

       SUBROUTINE ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,&
     &                    DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR, &
     &                    WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: B(LDB,*), DL(*), D(*), DU(*) 
         COMPLEX(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
      END SUBROUTINE ZGTSVX

      MODULE PROCEDURE DGTSVX1
      MODULE PROCEDURE ZGTSVX1

       END INTERFACE

      INTERFACE LA_GBSVX

       SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF,     &
     &                    LDAF, PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,&
     &                    FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE DGBSVX

       SUBROUTINE ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF,     &
     &                    LDAF, PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,&
     &                    FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE ZGBSVX

      MODULE PROCEDURE DGBSVX1
      MODULE PROCEDURE ZGBSVX1

       END INTERFACE

      INTERFACE LA_GESVX

       SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,  &
     &                    EQUED, R, C, B, LDB, X, LDX, RCOND, FERR,     &
     &                    BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE DGESVX

       SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,  &
     &                    EQUED, R, C, B, LDB, X, LDX, RCOND, FERR,     &
     &                    BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
      END SUBROUTINE ZGESVX

      MODULE PROCEDURE DGESVX1
      MODULE PROCEDURE ZGESVX1

      END INTERFACE

      INTERFACE LA_GETRF

       SUBROUTINE DGETRF( M, N, A, LDA, PIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT( OUT ) :: PIV( * )
         REAL(WP), INTENT( INOUT ) :: A( LDA, * )
      END SUBROUTINE DGETRF

       SUBROUTINE ZGETRF( M, N, A, LDA, PIV, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT( OUT ) :: PIV( * )
         COMPLEX(WP), INTENT( INOUT ) :: A( LDA, * )
      END SUBROUTINE ZGETRF


       END INTERFACE

      INTERFACE LA_GETRS

       SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: PIV(*)
         REAL(WP), INTENT(IN) :: A(LDA,*)
         REAL(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE DGETRS

       SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: PIV(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
      END SUBROUTINE ZGETRS

      MODULE PROCEDURE DGETRS1
      MODULE PROCEDURE ZGETRS1

       END INTERFACE

      INTERFACE LA_GETRI

       SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE DGETRI

       SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE ZGETRI

       END INTERFACE

      INTERFACE LA_GERFS

       SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B, LDB,&
     &                    X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: PIV(*)
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: X(LDX,*)
      END SUBROUTINE DGERFS

       SUBROUTINE ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B, LDB,&
     &                    X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: PIV(*)
         REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
      END SUBROUTINE ZGERFS

      MODULE PROCEDURE DGERFS1
      MODULE PROCEDURE ZGERFS1

       END INTERFACE

      INTERFACE LA_GEEQU

       SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,     &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, COLCND, ROWCND
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: C( * ), R( * )
      END SUBROUTINE DGEEQU

       SUBROUTINE ZGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,     &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, M, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: AMAX, COLCND, ROWCND
         COMPLEX(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: C( * ), R( * )
      END SUBROUTINE ZGEEQU

       END INTERFACE

      INTERFACE LA_LANGE

       FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLANGE
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, M, N
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION DLANGE

       FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: ZLANGE
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, M, N
         COMPLEX(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION ZLANGE

      MODULE PROCEDURE DLANGE1
      MODULE PROCEDURE ZLANGE1

       END INTERFACE

      INTERFACE LA_GECON

       SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DGECON

       SUBROUTINE ZGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, RWORK,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA, * )
         COMPLEX(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZGECON

       END INTERFACE

      INTERFACE LA_SYEV

       SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYEV

       END INTERFACE

      INTERFACE LA_HEEV

       SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  &
     &                   INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHEEV

       END INTERFACE

      INTERFACE LA_SYEVD

       SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
     &                    LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDA, LIWORK, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYEVD

       END INTERFACE

      INTERFACE LA_HEEVD

       SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
     &                    LRWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: LDA, LIWORK, LWORK, N, LRWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE ZHEEVD

       END INTERFACE

       INTERFACE LA_SYEVR

       SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,    &
     &                    IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: N, IL, IU, LDZ, LDA, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: M
         INTEGER, INTENT(OUT) :: ISUPPZ(*)
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
       END SUBROUTINE  DSYEVR    

      END INTERFACE
      
      INTERFACE LA_HEEVR
		   
       SUBROUTINE ZHEEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,    &
     &                    RWORK, LRWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: N, IL, IU, LDZ, LDA, LWORK, LIWORK, LRWORK
         INTEGER, INTENT(OUT) :: M
         INTEGER, INTENT(OUT) :: ISUPPZ(*)
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
        END SUBROUTINE ZHEEVR 


       END INTERFACE

      INTERFACE LA_SYEVX

       SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,     &
     &                    IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: IL, IU, LDA, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, M
         INTEGER, INTENT(OUT) :: IFAIL(*), IWORK(*)
         REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
      END SUBROUTINE DSYEVX

       END INTERFACE

      INTERFACE LA_HEEVX

       SUBROUTINE ZHEEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
     &                    ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK,     &
     &                    IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: IL, IU, LDA, LDZ, LWORK, N
         INTEGER, INTENT(OUT) :: INFO, M
         INTEGER, INTENT(OUT) :: IFAIL(*), IWORK(*)
         REAL(WP), INTENT(IN) :: ABSTOL, VL, VU
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
      END SUBROUTINE ZHEEVX

       END INTERFACE

      INTERFACE LA_SYGST

       SUBROUTINE DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: ITYPE, LDA, LDB, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: B(LDB,*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE DSYGST

       END INTERFACE

      INTERFACE LA_HEGST

       SUBROUTINE ZHEGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: ITYPE, LDA, LDB, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE ZHEGST

       END INTERFACE

      INTERFACE LA_SYGV

       SUBROUTINE DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
     &                   LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: ITYPE, LDA, LDB, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
      END SUBROUTINE DSYGV

       END INTERFACE

      INTERFACE LA_HEGV

       SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
     &                   LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: ITYPE, LDA, LDB, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: W(*), RWORK(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
       END SUBROUTINE ZHEGV

       END INTERFACE

        INTERFACE LA_SYGVX
	

       SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,  &
     &                    VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,   &
     &                    LWORK, IWORK, IFAIL, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
          INTEGER, INTENT(IN) :: ITYPE, N, IL, IU, LDZ, LDA, LDB, LWORK
          INTEGER, INTENT(OUT) :: M
          REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
          INTEGER, INTENT(OUT) ::  IWORK(*)
          INTEGER, INTENT(OUT) :: INFO
          REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
          REAL(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
          REAL(WP), INTENT(OUT) :: W(*)	
          INTEGER, INTENT(IN) :: IFAIL(*)
         END SUBROUTINE DSYGVX
   
        END INTERFACE
		
        INTERFACE LA_HEGVX
	
       SUBROUTINE ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,  &
     &                    VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,   &
     &                    LWORK, RWORK, IWORK, IFAIL, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, RANGE, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, IL, IU, LDZ, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: M
         REAL(WP), INTENT(IN) ::  ABSTOL, VL, VU
         INTEGER, INTENT(OUT) ::  IWORK(*)
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*), Z(LDZ,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         INTEGER, INTENT(IN) :: IFAIL(*)
       END SUBROUTINE ZHEGVX

       END INTERFACE


       INTERFACE LA_SYGVD


       SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,&
     &                    LWORK, IWORK, LIWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
         INTEGER, INTENT(IN) :: ITYPE, N, LDA, LDB, LWORK, LIWORK
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         REAL(WP), INTENT(OUT) :: W(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
        END SUBROUTINE DSYGVD
	
        END INTERFACE

        INTERFACE LA_HEGVD
	
       SUBROUTINE ZHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,&
     &                    LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          CHARACTER(LEN=1), INTENT(IN) :: JOBZ, UPLO
          INTEGER, INTENT(IN) :: ITYPE, N, LDA, LDB, LWORK, LIWORK, LRWORK
          INTEGER, INTENT(OUT) :: INFO, IWORK(*)
          COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
          REAL(WP), INTENT(OUT) :: W(*)
          COMPLEX(WP), INTENT(OUT) :: WORK(*)
          REAL(WP), INTENT(OUT) :: RWORK(*)
         END SUBROUTINE ZHEGVD
	
      END INTERFACE
            
      INTERFACE LA_SYTRD

       SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK,      &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         REAL(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE DSYTRD

       END INTERFACE

      INTERFACE LA_HETRD

       SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK,      &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: D(*), E(*)
         COMPLEX(WP), INTENT(OUT) :: TAU(*), WORK(LWORK)
      END SUBROUTINE ZHETRD

       END INTERFACE

      INTERFACE LA_ORGTR

       SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: TAU(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
         REAL(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE DORGTR

       END INTERFACE

      INTERFACE LA_UNGTR

       SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LWORK, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: TAU(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(LWORK)
      END SUBROUTINE ZUNGTR

       END INTERFACE

      INTERFACE LA_LANSY

      FUNCTION DLANSY( NORM, UPLO, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLANSY
         CHARACTER(LEN=1), INTENT(IN) :: NORM, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION DLANSY

      FUNCTION ZLANSY( NORM, UPLO, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: ZLANSY
         CHARACTER(LEN=1), INTENT(IN) :: NORM, UPLO
         INTEGER, INTENT(IN) :: LDA, N
         COMPLEX(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END FUNCTION ZLANSY

       END INTERFACE

      INTERFACE LA_POTRF

       SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE DPOTRF

       SUBROUTINE ZPOTRF( UPLO, N, A, LDA, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*)
      END SUBROUTINE ZPOTRF

       END INTERFACE

      INTERFACE LA_POCON

       SUBROUTINE DPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, IWORK,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND
         INTEGER, INTENT(OUT) :: IWORK( * )
         REAL(WP), INTENT(IN) :: A( LDA, * )
         REAL(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE DPOCON

       SUBROUTINE ZPOCON( UPLO, N, A, LDA, ANORM, RCOND, WORK, RWORK,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: ANORM
         REAL(WP), INTENT(OUT) :: RCOND, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA, * )
         COMPLEX(WP), INTENT(OUT) :: WORK( * )
      END SUBROUTINE ZPOCON

      END INTERFACE

      INTERFACE LA_ILAENV

      FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
         INTEGER :: ILAENV
         CHARACTER(LEN=*), INTENT(IN) :: NAME, OPTS
         INTEGER, INTENT(IN) :: ISPEC, N1, N2, N3, N4
      END FUNCTION ILAENV

      END INTERFACE

      INTERFACE LA_LAGGE

       SUBROUTINE DLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, KL, KU, LDA
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: ISEED(4)
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(OUT) :: A(LDA,*), WORK(*)
      END SUBROUTINE DLAGGE

       SUBROUTINE ZLAGGE( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, KL, KU, LDA
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: ISEED(4)
         REAL(WP), INTENT(IN) :: D(*)
         COMPLEX(WP), INTENT(OUT) :: A(LDA,*), WORK(*)
      END SUBROUTINE ZLAGGE

      END INTERFACE

      CONTAINS

      SUBROUTINE DGESV1( N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         INTERFACE
           SUBROUTINE DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: PIV(*)
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           END SUBROUTINE DGESV
         END INTERFACE
         CALL DGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
      END SUBROUTINE DGESV1

      SUBROUTINE ZGESV1( N, NRHS, A, LDA, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         INTERFACE
           SUBROUTINE ZGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: PIV(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           END SUBROUTINE ZGESV
         END INTERFACE
         CALL ZGESV( N, NRHS, A, LDA, PIV, B, LDB, INFO )
      END SUBROUTINE ZGESV1



      SUBROUTINE DGESVX1( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,  &
     &                    EQUED, R, C, B, LDB, X, LDX, RCOND, FERR,     &
     &                    BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF,   &
     &                        PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,  &
     &                        FERR, BERR, WORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             INTEGER, INTENT(INOUT) :: PIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: R(*), C(*)
             REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
           END SUBROUTINE DGESVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,      &
     &                EQUED, R, C, B, LDB, X, LDX, RCOND, LFERR, LBERR, &
     &                WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DGESVX1

      SUBROUTINE ZGESVX1( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,  &
     &                    EQUED, R, C, B, LDB, X, LDX, RCOND, FERR,     &
     &                    BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF,   &
     &                        PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,  &
     &                        FERR, BERR, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: PIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: R(*), C(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*),        &
     &                                     B(LDB,*)
           END SUBROUTINE ZGESVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,      &
     &                EQUED, R, C, B, LDB, X, LDX, RCOND, LFERR, LBERR, &
     &                WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZGESVX1


      SUBROUTINE DPOSV1( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         INTERFACE
           SUBROUTINE DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           END SUBROUTINE DPOSV
         END INTERFACE
         CALL DPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      END SUBROUTINE DPOSV1

      SUBROUTINE ZPOSV1( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         INTERFACE
           SUBROUTINE ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB, LDA
             INTEGER, INTENT(OUT) :: INFO
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           END SUBROUTINE ZPOSV
         END INTERFACE
         CALL ZPOSV( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      END SUBROUTINE ZPOSV1


      FUNCTION DLANGE1( NORM, M, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: DLANGE1
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, M, N
         REAL(WP), INTENT(IN) :: A( * )
         REAL(WP), INTENT(OUT) :: WORK( * )
         INTERFACE
           FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
             USE LA_PRECISION, ONLY: WP => DP
             REAL(WP) :: DLANGE
             CHARACTER(LEN=1), INTENT(IN) :: NORM
             INTEGER, INTENT(IN) :: LDA, M, N
             REAL(WP), INTENT(IN) :: A( LDA, * )
             REAL(WP), INTENT(OUT) :: WORK( * )
           END FUNCTION DLANGE
         END INTERFACE
        DLANGE1 = DLANGE( NORM, M, N, A, LDA, WORK )
      END FUNCTION DLANGE1

      FUNCTION ZLANGE1( NORM, M, N, A, LDA, WORK )
         USE LA_PRECISION, ONLY: WP => DP
         REAL(WP) :: ZLANGE1
         CHARACTER(LEN=1), INTENT(IN) :: NORM
         INTEGER, INTENT(IN) :: LDA, M, N
         COMPLEX(WP), INTENT(IN) :: A( * )
         REAL(WP), INTENT(OUT) :: WORK( * )
         INTERFACE
           FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
             USE LA_PRECISION, ONLY: WP => DP
             REAL(WP) :: ZLANGE
             CHARACTER(LEN=1), INTENT(IN) :: NORM
             INTEGER, INTENT(IN) :: LDA, M, N
             COMPLEX(WP), INTENT(IN) :: A( LDA, * )
             REAL(WP), INTENT(OUT) :: WORK( * )
           END FUNCTION ZLANGE
         END INTERFACE
        ZLANGE1 = ZLANGE( NORM, M, N, A, LDA, WORK )
      END FUNCTION ZLANGE1


      SUBROUTINE DGBSV1( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(*)
         INTERFACE
           SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB,    &
     &                       INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: PIV(*)
             REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
           END SUBROUTINE DGBSV
         END INTERFACE
         CALL DGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
      END SUBROUTINE DGBSV1

      SUBROUTINE ZGBSV1( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: PIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(*)
         INTERFACE
           SUBROUTINE ZGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB,    &
     &                       INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: PIV(*)
             COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
           END SUBROUTINE ZGBSV
         END INTERFACE
         CALL ZGBSV( N, KL, KU, NRHS, AB, LDAB, PIV, B, LDB, INFO )
      END SUBROUTINE ZGBSV1


      SUBROUTINE DGBSVX1( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF,     &
     &                    LDAF, PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,&
     &                    FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE DGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF, &
     &                        LDAF, PIV, EQUED, R, C, B, LDB, X, LDX,   &
     &                        RCOND, FERR, BERR, WORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
             INTEGER, INTENT(OUT) :: INFO
!
             INTEGER, INTENT(OUT) :: IWORK(*)
             INTEGER, INTENT(INOUT) :: PIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: R(*), C(*)
             REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
           END SUBROUTINE DGBSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF, LDAF,   &
     &                PIV, EQUED, R, C, B, LDB, X, LDX, RCOND, LFERR,   &
     &                LBERR, WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DGBSVX1

      SUBROUTINE ZGBSVX1( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF,     &
     &                    LDAF, PIV, EQUED, R, C, B, LDB, X, LDX, RCOND,&
     &                    FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: PIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: R(*), C(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF, &
     &                        LDAF, PIV, EQUED, R, C, B, LDB, X, LDX,   &
     &                        RCOND, FERR, BERR, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, KL, KU
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: PIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: R(*), C(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
           END SUBROUTINE ZGBSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZGBSVX( FACT, TRANS, N, KL, KU, NRHS, A, LDA, AF, LDAF,   &
     &                PIV, EQUED, R, C, B, LDB, X, LDX, RCOND, LFERR,   &
     &                LBERR, WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZGBSVX1


      SUBROUTINE DGTSV1( N, NRHS, DL, D, DU, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(*)
         INTERFACE
           SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(LDB,*)
           END SUBROUTINE DGTSV
         END INTERFACE
         CALL DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
      END SUBROUTINE DGTSV1

      SUBROUTINE ZGTSV1( N, NRHS, DL, D, DU, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(*)
         INTERFACE
           SUBROUTINE ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             COMPLEX(WP), INTENT(INOUT) :: DL(*), D(*), DU(*), B(LDB,*)
           END SUBROUTINE ZGTSV
         END INTERFACE
         CALL ZGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
      END SUBROUTINE ZGTSV1


       SUBROUTINE DGTSVX1( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF,    &
     &                     DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, &
     &                     BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, X(*), WORK(*)
         REAL(WP), INTENT(IN) :: B(*), DL(*), D(*), DU(*) 
         REAL(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
         INTERFACE
           SUBROUTINE DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, &
     &                        DUF, DU2, IPIV, B, LDB, X, LDX, RCOND,    &
     &                        FERR, BERR, WORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), X(LDX,*), WORK(*)
             REAL(WP), INTENT(IN) :: B(LDB,*), DL(*), D(*), DU(*) 
             REAL(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
           END SUBROUTINE DGTSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,    &
     &                DU2, IPIV, B, LDB, X, LDX, RCOND, LFERR, LBERR,   &
     &                WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DGTSVX1

       SUBROUTINE ZGTSVX1( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF,    &
     &                     DUF, DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, &
     &                     BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: B(*), DL(*), D(*), DU(*) 
         COMPLEX(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
         INTERFACE
           SUBROUTINE ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, &
     &                        DUF, DU2, IPIV, B, LDB, X, LDX, RCOND,    &
     &                        FERR, BERR, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS, FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             COMPLEX(WP), INTENT(IN) :: B(LDB,*), DL(*), D(*), DU(*) 
             COMPLEX(WP), INTENT(INOUT) :: DF(*), DLF(*), DU2(*), DUF(*)
           END SUBROUTINE ZGTSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,    &
     &                DU2, IPIV, B, LDB, X, LDX, RCOND, LFERR, LBERR,   &
     &                WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZGTSVX1


       SUBROUTINE DPOSVX1( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED,&
     &                     S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,  &
     &                     IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,    &
     &                        EQUED, S, B, LDB, X, LDX, RCOND, FERR,    &
     &                        BERR, WORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: S(*)
             REAL(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
           END SUBROUTINE DPOSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S,  &
     &                B, LDB, X, LDX, RCOND, LFERR, LBERR, WORK, IWORK, &
     &                INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DPOSVX1

       SUBROUTINE ZPOSVX1( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED,&
     &                     S, B, LDB, X, LDX, RCOND, FERR, BERR, WORK,  &
     &                     RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(*)
         INTERFACE
           SUBROUTINE ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,    &
     &                        EQUED, S, B, LDB, X, LDX, RCOND, FERR,    &
     &                        BERR, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: S(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
           END SUBROUTINE ZPOSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZPOSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED, S,  &
     &                B, LDB, X, LDX, RCOND, LFERR, LBERR, WORK, RWORK, &
     &                INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZPOSVX1


       SUBROUTINE DPPSV1( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AP(*), B(*)
         INTERFACE
           SUBROUTINE DPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
           END SUBROUTINE DPPSV
         END INTERFACE
         CALL DPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE DPPSV1

       SUBROUTINE ZPPSV1( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(*)
         INTERFACE
           SUBROUTINE ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
           END SUBROUTINE ZPPSV
         END INTERFACE
         CALL ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE ZPPSV1


       SUBROUTINE DPPSVX1( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,   &
     &                     LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, &
     &                     INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: AP(*), AFP(*), B(*), S(*)
         INTERFACE
           SUBROUTINE DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,&
     &                        LDB, X, LDX, RCOND, FERR, BERR, WORK,     &
     &                        IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: AP(*), AFP(*), B(LDB,*), S(*)
           END SUBROUTINE DPPSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X,&
     &                LDX, RCOND, LFERR, LBERR, WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DPPSVX1

       SUBROUTINE ZPPSVX1( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,   &
     &                     LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, &
     &                     INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), AFP(*), B(*)
         INTERFACE
           SUBROUTINE ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B,&
     &                        LDB, X, LDX, RCOND, FERR, BERR, WORK,     &
     &                        RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: S(*)
             COMPLEX(WP), INTENT(INOUT) :: AP(*), AFP(*), B(LDB,*)
           END SUBROUTINE ZPPSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZPPSVX( FACT, UPLO, N, NRHS, AP, AFP, EQUED, S, B, LDB, X,&
     &                LDX, RCOND, LFERR, LBERR, WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZPPSVX1


       SUBROUTINE DPBSV1( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(*)
         INTERFACE
           SUBROUTINE DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB,       &
     &                       INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
           END SUBROUTINE DPBSV
         END INTERFACE
         CALL DPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
      END SUBROUTINE DPBSV1

       SUBROUTINE ZPBSV1( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(*)
         INTERFACE
           SUBROUTINE ZPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB,       &
     &                       INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB, KD, LDAB
             INTEGER, INTENT(OUT) :: INFO
             COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), B(LDB,*)
           END SUBROUTINE ZPBSV
         END INTERFACE
         CALL ZPBSV( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
      END SUBROUTINE ZPBSV1


       SUBROUTINE DPBSVX1( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB,      &
     &                     LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR,&
     &                     BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: X(*), WORK(*), FERR, BERR, RCOND
         REAL(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*), B(*),     &
     &                              S(*)
         INTERFACE
           SUBROUTINE DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB,   &
     &                        LDAFB, EQUED, S, B, LDB, X, LDX, RCOND,   &
     &                        FERR, BERR, WORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
             INTEGER, INTENT(OUT) :: INFO, IWORK(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*), FERR(*),       &
     &                                BERR(*), RCOND
             REAL(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*),       &
     &                                  B(LDB,*), S(*)
           END SUBROUTINE DPBSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB,    &
     &                EQUED, S, B, LDB, X, LDX, RCOND, LFERR, LBERR,    &
     &                WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DPBSVX1

       SUBROUTINE ZPBSVX1( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB,      &
     &                     LDAFB, EQUED, S, B, LDB, X, LDX, RCOND, FERR,&
     &                     BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
         CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
         INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: FERR, BERR, RCOND, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(INOUT) :: S(*)
         COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*), B(*)
         INTERFACE
           SUBROUTINE ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB,   &
     &                        LDAFB, EQUED, S, B, LDB, X, LDX, RCOND,   &
     &                        FERR, BERR, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT, UPLO
             CHARACTER(LEN=1), INTENT(INOUT) :: EQUED
             INTEGER, INTENT(IN) :: LDAB, LDAFB, LDB, LDX, NRHS, N, KD
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RCOND, RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(INOUT) :: S(*)
             COMPLEX(WP), INTENT(INOUT) :: AB(LDAB,*), AFB(LDAFB,*),    &
     &                                     B(LDB,*)
           END SUBROUTINE ZPBSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZPBSVX( FACT, UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB,    &
     &                EQUED, S, B, LDB, X, LDX, RCOND, LFERR, LBERR,    &
     &                WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZPBSVX1


       SUBROUTINE DPTSV1( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*)
         REAL(WP), INTENT(INOUT) :: E(*), B(*)
         INTERFACE
           SUBROUTINE DPTSV( N, NRHS, D, E, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: D(*)
             REAL(WP), INTENT(INOUT) :: E(*), B(LDB,*)
           END SUBROUTINE DPTSV
         END INTERFACE
         CALL DPTSV( N, NRHS, D, E, B, LDB, INFO )
      END SUBROUTINE DPTSV1

       SUBROUTINE ZPTSV1( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: D(*)
         COMPLEX(WP), INTENT(INOUT) :: E(*), B(*)
         INTERFACE
           SUBROUTINE ZPTSV( N, NRHS, D, E, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: D(*)
             COMPLEX(WP), INTENT(INOUT) :: E(*), B(LDB,*)
           END SUBROUTINE ZPTSV
         END INTERFACE
         CALL ZPTSV( N, NRHS, D, E, B, LDB, INFO )
      END SUBROUTINE ZPTSV1


       SUBROUTINE DPTSVX1( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, &
     &                     RCOND, FERR, BERR, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(IN) :: E(*), B(*)
         REAL(WP), INTENT(INOUT) :: DF(*)
         REAL(WP), INTENT(INOUT) :: EF(*)
         INTERFACE
           SUBROUTINE DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X,   &
     &                        LDX, RCOND, FERR, BERR, WORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(IN) :: D(*)
             REAL(WP), INTENT(IN) :: E(*), B(LDB,*)
             REAL(WP), INTENT(INOUT) :: DF(*)
             REAL(WP), INTENT(INOUT) :: EF(*)
           END SUBROUTINE DPTSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,      &
     &                RCOND, LFERR, LBERR, WORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DPTSVX1

       SUBROUTINE ZPTSVX1( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX, &
     &                     RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(IN) :: D(*)
         COMPLEX(WP), INTENT(IN) :: E(*), B(*)
         REAL(WP), INTENT(INOUT) :: DF(*)
         COMPLEX(WP), INTENT(INOUT) :: EF(*)
         INTERFACE
           SUBROUTINE ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X,   &
     &                        LDX, RCOND, FERR, BERR, WORK, RWORK,      &
     &                        INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(IN) :: D(*)
             COMPLEX(WP), INTENT(IN) :: E(*), B(LDB,*)
             REAL(WP), INTENT(INOUT) :: DF(*)
             COMPLEX(WP), INTENT(INOUT) :: EF(*)
           END SUBROUTINE ZPTSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,      &
     &                RCOND, LFERR, LBERR, WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZPTSVX1


       SUBROUTINE DSYSV1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,    &
     &                    LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
     &                       LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE DSYSV
         END INTERFACE
         CALL DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK,  &
     &               INFO )
      END SUBROUTINE DSYSV1

       SUBROUTINE ZSYSV1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,    &
     &                    LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
     &                       LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             COMPLEX(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE ZSYSV
         END INTERFACE
         CALL ZSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK,  &
     &               INFO )
      END SUBROUTINE ZSYSV1


       SUBROUTINE ZHESV1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,    &
     &                    LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, &
     &                       LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             COMPLEX(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE ZHESV
         END INTERFACE
         CALL ZHESV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK,  &
     &               INFO )
      END SUBROUTINE ZHESV1


       SUBROUTINE DSYSVX1( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, &
     &                     B, LDB, X, LDX, RCOND, FERR, BERR, WORK,     &
     &                     LWORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(IN) :: A(LDA,*), B(*)
         REAL(WP), INTENT(INOUT) :: AF(LDAF,*)
         INTERFACE
           SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,    &
     &                        IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,  &
     &                        WORK, LWORK, IWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(INOUT) :: AF(LDAF,*)
           END SUBROUTINE DSYSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,   &
     &                LDB, X, LDX, RCOND, LFERR, LBERR, WORK, LWORK,    &
     &                IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DSYSVX1

       SUBROUTINE ZSYSVX1( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, &
     &                     B, LDB, X, LDX, RCOND, FERR, BERR, WORK,     &
     &                     LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
         INTERFACE
           SUBROUTINE ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,    &
     &                        IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,  &
     &                        WORK, LWORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
             COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
           END SUBROUTINE ZSYSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,   &
     &                LDB, X, LDX, RCOND, LFERR, LBERR, WORK, LWORK,    &
     &                RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZSYSVX1


       SUBROUTINE ZHESVX1( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, &
     &                     B, LDB, X, LDX, RCOND, FERR, BERR, WORK,     &
     &                     LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
         INTERFACE
           SUBROUTINE ZHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,    &
     &                        IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,  &
     &                        WORK, LWORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N, LWORK
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
             COMPLEX(WP), INTENT(INOUT) :: AF(LDAF,*)
      END SUBROUTINE ZHESVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZHESVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,   &
     &                LDB, X, LDX, RCOND, LFERR, LBERR, WORK, LWORK,    &
     &                RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZHESVX1


       SUBROUTINE DSPSV1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: AP(*), B(*)
         INTERFACE
           SUBROUTINE DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             REAL(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
           END SUBROUTINE DSPSV
         END INTERFACE
         CALL DSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE DSPSV1

       SUBROUTINE ZSPSV1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(*)
         INTERFACE
           SUBROUTINE ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
           END SUBROUTINE ZSPSV
         END INTERFACE
         CALL ZSPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE ZSPSV1


       SUBROUTINE ZHPSV1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: NRHS, N, LDB
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AP(*), B(*)
         INTERFACE
           SUBROUTINE ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO
             INTEGER, INTENT(IN) :: NRHS, N, LDB
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(IN) :: IPIV(*)
             COMPLEX(WP), INTENT(INOUT) :: AP(*), B(LDB,*)
           END SUBROUTINE ZHPSV
         END INTERFACE
         CALL ZHPSV( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE ZHPSV1


       SUBROUTINE DSPSVX1( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, &
     &                     LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IWORK(*)
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR
         REAL(WP), INTENT(OUT) :: X(*), WORK(*)
         REAL(WP), INTENT(IN) :: A(*), B(*)
         REAL(WP), INTENT(INOUT) :: AF(*)
         INTERFACE
           SUBROUTINE DSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, &
     &                        X, LDX, RCOND, FERR, BERR, WORK, IWORK,   &
     &                        INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(OUT) :: IWORK(*)
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
             REAL(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             REAL(WP), INTENT(IN) :: A(*), B(LDB,*)
             REAL(WP), INTENT(INOUT) :: AF(*)
           END SUBROUTINE DSPSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL DSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, LDX, &
     &                RCOND, LFERR, LBERR, WORK, IWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE DSPSVX1

       SUBROUTINE ZSPSVX1( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, &
     &                     LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: AF(*)
         INTERFACE
           SUBROUTINE ZSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, &
     &                        X, LDX, RCOND, FERR, BERR, WORK, RWORK,   &
     &                        INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             COMPLEX(WP), INTENT(IN) :: A(*), B(LDB,*)
             COMPLEX(WP), INTENT(INOUT) :: AF(*)
           END SUBROUTINE ZSPSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZSPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, LDX, &
     &                RCOND, LFERR, LBERR, WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZSPSVX1


       SUBROUTINE ZHPSVX1( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, &
     &                     LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
         INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: RCOND
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: X(*), WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: AF(*)
         INTERFACE
           SUBROUTINE ZHPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, &
     &                        X, LDX, RCOND, FERR, BERR, WORK, RWORK,   &
     &                        INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: UPLO, FACT
             INTEGER, INTENT(IN) :: LDB, LDX, NRHS, N
             INTEGER, INTENT(OUT) :: INFO
             INTEGER, INTENT(INOUT) :: IPIV(*)
             REAL(WP), INTENT(OUT) :: RCOND
             REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: X(LDX,*), WORK(*)
             COMPLEX(WP), INTENT(IN) :: A(*), B(LDB,*)
             COMPLEX(WP), INTENT(INOUT) :: AF(*)
           END SUBROUTINE ZHPSVX
         END INTERFACE
         REAL(WP) :: LFERR(1), LBERR(1)
         CALL ZHPSVX( FACT, UPLO, N, NRHS, A, AF, IPIV, B, LDB, X, LDX, &
     &                RCOND, LFERR, LBERR, WORK, RWORK, INFO )
         FERR = LFERR(1); BERR = LBERR(1)
      END SUBROUTINE ZHPSVX1


       SUBROUTINE DGELS1( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK,      &
     &                    LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK,   &
     &                       LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE DGELS
         END INTERFACE
         CALL DGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,    &
     &               INFO )
      END SUBROUTINE DGELS1

       SUBROUTINE ZGELS1( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK,      &
     &                    LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK,   &
     &                       LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             CHARACTER(LEN=1), INTENT(IN) :: TRANS
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             COMPLEX(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE ZGELS
         END INTERFACE
         CALL ZGELS( TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK,    &
     &               INFO )
      END SUBROUTINE ZGELS1


       SUBROUTINE DGELSY1( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,     &
     &                     RANK, WORK, LWORK, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
          INTEGER, INTENT(OUT) :: INFO, RANK
          INTEGER, INTENT(INOUT) :: JPVT(*)
          REAL(WP), INTENT(IN) :: RCOND
          REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
          REAL(WP), INTENT(OUT) ::  WORK(*)
          INTERFACE
           SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,  &
     &                        RANK, WORK, LWORK, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
          INTEGER, INTENT(OUT) :: INFO, RANK
          INTEGER, INTENT(INOUT) :: JPVT(*)
          REAL(WP), INTENT(IN) :: RCOND
          REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
          REAL(WP), INTENT(OUT) ::  WORK(*)
         END SUBROUTINE DGELSY
        END INTERFACE
         CALL DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,    &
     &                WORK, LWORK, INFO )
       END SUBROUTINE DGELSY1
       SUBROUTINE ZGELSY1( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,     &
     &                     RANK, WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(OUT) ::  WORK(*)
         REAL(WP) :: RWORK(*)
         INTERFACE
           SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,  &
     &                        RANK, WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: M, N, NRHS, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
         COMPLEX(WP), INTENT(OUT) ::  WORK(*)
         REAL(WP) :: RWORK(*)
        END SUBROUTINE ZGELSY
       END INTERFACE
         CALL ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,    &
     &                WORK, LWORK, RWORK, INFO )
       END SUBROUTINE ZGELSY1

       SUBROUTINE DGELSD1( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP =>  DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND 
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: S(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTEGER :: IWORK(*)  
         INTERFACE
           SUBROUTINE DGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND,     &
     &                        RANK, WORK, LWORK, IWORK, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
           INTEGER, INTENT(OUT) :: INFO, RANK
           REAL(WP), INTENT(IN) :: RCOND
           REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           REAL(WP), INTENT(OUT) :: S(*)
           REAL(WP), INTENT(OUT) :: WORK(*)
           INTEGER :: IWORK(*)
         END SUBROUTINE DGELSD
       END INTERFACE
         CALL DGELSD ( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK,&
     &                 LWORK, IWORK, INFO )
       END SUBROUTINE DGELSD1
       SUBROUTINE ZGELSD1( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, RWORK, IWORK, INFO )
       USE LA_PRECISION, ONLY: WP =>  DP
       INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
       INTEGER, INTENT(OUT) :: INFO, RANK
       REAL(WP), INTENT(IN) :: RCOND
       COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
       REAL(WP), INTENT(OUT) :: S(*)
       COMPLEX(WP), INTENT(OUT) :: WORK(*)
       INTEGER :: IWORK(*)
       REAL(WP) :: RWORK(*)
       INTERFACE
           SUBROUTINE ZGELSD( M, N, NRHS, A, LDA, B, LDB, S, RCOND,     &
     &                        RANK, WORK, LWORK, RWORK, IWORK, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
           INTEGER, INTENT(OUT) :: INFO, RANK
           REAL(WP), INTENT(IN) :: RCOND
           COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
           REAL(WP), INTENT(OUT) :: S(*)
           COMPLEX(WP), INTENT(OUT) :: WORK(*)
           INTEGER :: IWORK(*)
           REAL(WP) :: RWORK(*)
         END SUBROUTINE ZGELSD
       END INTERFACE
         CALL ZGELSD ( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK,&
     &                 LWORK, RWORK, IWORK, INFO )
      END SUBROUTINE ZGELSD1


       SUBROUTINE DGELSX1( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,     &
     &                     RANK, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE DGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,  &
     &                        RANK, WORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
             INTEGER, INTENT(OUT) :: INFO, RANK
             INTEGER, INTENT(INOUT) :: JPVT(*)
             REAL(WP), INTENT(IN) :: RCOND
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE DGELSX
         END INTERFACE
         CALL DGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,    &
     &                WORK, INFO )
      END SUBROUTINE DGELSX1

       SUBROUTINE ZGELSX1( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,     &
     &                     RANK, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
         INTEGER, INTENT(OUT) :: INFO, RANK
         INTEGER, INTENT(INOUT) :: JPVT(*)
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND,  &
     &                        RANK, WORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB
             INTEGER, INTENT(OUT) :: INFO, RANK
             INTEGER, INTENT(INOUT) :: JPVT(*)
             REAL(WP), INTENT(IN) :: RCOND
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE ZGELSX
         END INTERFACE
         CALL ZGELSX( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK,    &
     &                WORK, RWORK, INFO )
      END SUBROUTINE ZGELSX1


       SUBROUTINE DGELSS1( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND
         REAL(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: S(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND,     &
     &                        RANK, WORK, LWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO, RANK
             REAL(WP), INTENT(IN) :: RCOND
             REAL(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: S(*)
             REAL(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE DGELSS
         END INTERFACE
         CALL DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, &
     &                LWORK, INFO )
      END SUBROUTINE DGELSS1

       SUBROUTINE ZGELSS1( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,  &
     &                     WORK, LWORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
         INTEGER, INTENT(OUT) :: INFO, RANK
         REAL(WP), INTENT(IN) :: RCOND
         COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(*)
         REAL(WP), INTENT(OUT) :: S(*), RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
           SUBROUTINE ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND,     &
     &                        RANK, WORK, LWORK, RWORK, INFO )
             USE LA_PRECISION, ONLY: WP => DP
             INTEGER, INTENT(IN) :: NRHS, M, N, LDA, LDB, LWORK
             INTEGER, INTENT(OUT) :: INFO, RANK
             REAL(WP), INTENT(IN) :: RCOND
             COMPLEX(WP), INTENT(INOUT) :: A(LDA,*), B(LDB,*)
             REAL(WP), INTENT(OUT) :: S(*), RWORK(*)
             COMPLEX(WP), INTENT(OUT) :: WORK(*)
           END SUBROUTINE ZGELSS
         END INTERFACE
         CALL ZGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, &
     &                LWORK, RWORK, INFO )
      END SUBROUTINE ZGELSS1

       SUBROUTINE DGETRS1( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          CHARACTER(LEN=1), INTENT(IN) :: TRANS
          INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
          INTEGER, INTENT(OUT) :: INFO
          INTEGER, INTENT(IN) :: PIV(*)
          REAL(WP), INTENT(IN) :: A(LDA,*)
          REAL(WP), INTENT(INOUT) :: B(*)
          INTERFACE
             SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB,    &
     &                          INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: PIV(*)
               REAL(WP), INTENT(IN) :: A(LDA,*)
               REAL(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE DGETRS
          END INTERFACE
          CALL DGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
       END SUBROUTINE DGETRS1
       SUBROUTINE ZGETRS1( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
          USE LA_PRECISION, ONLY: WP => DP
          CHARACTER(LEN=1), INTENT(IN) :: TRANS
          INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
          INTEGER, INTENT(OUT) :: INFO
          INTEGER, INTENT(IN) :: PIV(*)
          COMPLEX(WP), INTENT(IN) :: A(LDA,*)
          COMPLEX(WP), INTENT(INOUT) :: B(*)
          INTERFACE
             SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB,    &
     &                          INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDA, LDB, NRHS, N
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: PIV(*)
               COMPLEX(WP), INTENT(IN) :: A(LDA,*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZGETRS
          END INTERFACE
          CALL ZGETRS( TRANS, N, NRHS, A, LDA, PIV, B, LDB, INFO )
       END SUBROUTINE ZGETRS1

       SUBROUTINE DGERFS1( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B,    &
     &                     LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
           USE LA_PRECISION, ONLY: WP => DP
           CHARACTER(LEN=1), INTENT(IN) :: TRANS
           INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
           INTEGER, INTENT(OUT) :: INFO
           INTEGER, INTENT(IN) :: PIV(*)
           INTEGER, INTENT(OUT) :: IWORK(*)
           REAL(WP), INTENT(OUT) :: FERR, BERR
           REAL(WP), INTENT(OUT) :: WORK(*)
           REAL(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(*)
           REAL(WP), INTENT(INOUT) :: X(*)
           INTERFACE
              SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, &
     &                           B, LDB, X, LDX, FERR, BERR, WORK,      &
     &                           IWORK, INFO )
                 USE LA_PRECISION, ONLY: WP => DP
                 CHARACTER(LEN=1), INTENT(IN) :: TRANS
                 INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
                 INTEGER, INTENT(OUT) :: INFO
                 INTEGER, INTENT(IN) :: PIV(*)
                 INTEGER, INTENT(OUT) :: IWORK(*)
                 REAL(WP), INTENT(OUT) :: FERR(*), BERR(*)
                 REAL(WP), INTENT(OUT) :: WORK(*)
                 REAL(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
                 REAL(WP), INTENT(INOUT) :: X(LDX,*)
              END SUBROUTINE DGERFS
           END INTERFACE
           REAL(WP) FERR1(1), BERR1(1)
           CALL DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B, LDB,  &
     &                  X, LDX, FERR1, BERR1, WORK, IWORK, INFO )
           FERR = FERR1(1); BERR = BERR1(1)
        END SUBROUTINE DGERFS1

       SUBROUTINE ZGERFS1( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B,    &
     &                     LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: PIV(*)
         REAL(WP), INTENT(OUT) :: FERR, BERR, RWORK(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         INTERFACE
             SUBROUTINE ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV,  &
     &                          B, LDB, X, LDX, FERR, BERR, WORK, RWORK,&
     &                          INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, NRHS, N
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: PIV(*)
               REAL(WP), INTENT(OUT) :: FERR(*), BERR(*), RWORK(*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
               COMPLEX(WP), INTENT(IN) :: A(LDA,*), AF(LDAF,*), B(LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
            END SUBROUTINE ZGERFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
            CALL ZGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, PIV, B, LDB, &
     &                   X, LDX, FERR1, BERR1, WORK, RWORK, INFO )
           FERR = FERR1(1); BERR = BERR1(1)
        END SUBROUTINE ZGERFS1


      SUBROUTINE DGBTRS1( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B,    &
     &                    LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         REAL(WP), INTENT(INOUT) :: AB( LDAB,*), B(*)
         INTERFACE
            SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV,  &
     &                         B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(INOUT) :: IPIV(*)
               REAL(WP), INTENT(INOUT) :: AB( LDAB,*), B(LDB,*)
            END SUBROUTINE DGBTRS
         END INTERFACE
         CALL DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,   &
     &                INFO )
      END SUBROUTINE DGBTRS1
      SUBROUTINE ZGBTRS1( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B,    &
     &                    LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(INOUT) :: IPIV(*)
         COMPLEX(WP), INTENT(INOUT) :: AB( LDAB,*), B(*)
         INTERFACE
            SUBROUTINE ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV,  &
     &                         B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: KL, KU, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(INOUT) :: IPIV(*)
               COMPLEX(WP), INTENT(INOUT) :: AB( LDAB,*), B(LDB,*)
            END SUBROUTINE ZGBTRS
         END INTERFACE
         CALL ZGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,   &
     &                INFO )
      END SUBROUTINE ZGBTRS1

      SUBROUTINE DGBRFS1( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, &
     &                    IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         INTEGER, INTENT(OUT) :: IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*), B(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         INTERFACE
            SUBROUTINE DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,   &
     &                         LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, &
     &                         WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, &
     &                                NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               INTEGER, INTENT(OUT) :: IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*),      &
     &                                 B( LDB,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
            END SUBROUTINE DGBRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL DGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,     &
     &                IPIV, B, LDB, X, LDX, FERR1, BERR1, WORK, IWORK,  &
     &                INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DGBRFS1

      SUBROUTINE ZGBRFS1( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, &
     &                    IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*), B(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         INTERFACE
            SUBROUTINE ZGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB,   &
     &                         LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, &
     &                         WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: KL, KU, LDAB, LDAFB, LDB, LDX, N, &
     &                                NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AB( LDAB,*), AFB( LDAFB,*),   &
     &                                    B( LDB,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
            END SUBROUTINE ZGBRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL ZGBRFS( TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB,     &
     &                IPIV, B, LDB, X, LDX, FERR1, BERR1, WORK, RWORK,  &
     &                INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZGBRFS1

      SUBROUTINE DGTTRS1( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IPIV(*)
         REAL(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, &
     &                         LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(OUT) :: IPIV(*)
               REAL(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
               REAL(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE DGTTRS
         END INTERFACE
         CALL DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,     &
     &                INFO )
      END SUBROUTINE DGTTRS1
      SUBROUTINE ZGTTRS1( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(OUT) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, &
     &                         LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(OUT) :: IPIV(*)
               COMPLEX(WP), INTENT(IN) :: D(*), DL(*), DU(*), DU2(*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZGTTRS
         END INTERFACE
         CALL ZGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB,     &
     &                INFO )
      END SUBROUTINE ZGTTRS1

      SUBROUTINE DGTRFS1( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, &
     &                    IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: B(*), D(*), DF(*), DL(*), DLF(*),      &
     &                           DU(*), DU2(*), DUF(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, &
     &                         DU2, IPIV, B, LDB, X, LDX, FERR, BERR,   &
     &                         WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: B( LDB,*), D(*), DF(*), DL(*),   &
     &                                 DLF(*), DU(*), DU2(*), DUF(*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DGTRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL DGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,     &
     &                IPIV, B, LDB, X, LDX, FERR1, BERR1, WORK, IWORK,  &
     &                INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DGTRFS1

      SUBROUTINE ZGTRFS1( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2, &
     &                    IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK,&
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: TRANS
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: B(*), D(*), DF(*), DL(*), DLF(*),   &
     &                              DU(*), DU2(*), DUF(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, &
     &                         DU2, IPIV, B, LDB, X, LDX, FERR, BERR,   &
     &                         WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: TRANS
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS, IPIV(*)
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: B( LDB,*), D(*), DF(*), DL(*),&
     &                                    DLF(*), DU(*), DU2(*), DUF(*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZGTRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL ZGTRFS( TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF, DU2,     &
     &                IPIV, B, LDB, X, LDX, FERR1, BERR1, WORK, RWORK,  &
     &                INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZGTRFS1


      SUBROUTINE DPOTRS1( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A( LDA,*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: A( LDA,*)
               REAL(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE DPOTRS
         END INTERFACE
         CALL DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      END SUBROUTINE DPOTRS1

      SUBROUTINE ZPOTRS1( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: A( LDA,*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZPOTRS
         END INTERFACE
         CALL ZPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      END SUBROUTINE ZPOTRS1


      SUBROUTINE DPORFS1( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,   &
     &                    LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: A(*), AF( LDAF,*), B( LDB,*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, &
     &                         X, LDX, FERR, BERR, WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: A( LDA,*), AF( LDAF,*),          &
     &                                 B( LDB,*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DPORFS
         END INTERFACE
         REAL(WP) :: BERR1(1), FERR1(1)
         CALL DPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX,  &
     &                FERR1, BERR1, WORK, IWORK, INFO )
         BERR = BERR1(1); FERR = FERR1(1)
      END SUBROUTINE DPORFS1

      SUBROUTINE ZPORFS1( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X,   &
     &                    LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A(*), AF( LDAF,*), B( LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, &
     &                         X, LDX, FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: A( LDA,*), AF( LDAF,*),       &
     &                                    B( LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZPORFS
         END INTERFACE
         REAL(WP) :: BERR1(1), FERR1(1)
         CALL ZPORFS( UPLO, N, NRHS, A, LDA, AF, LDAF, B, LDB, X, LDX,  &
     &                FERR1, BERR1, WORK, RWORK, INFO )
         BERR = BERR1(1); FERR = FERR1(1)
      END SUBROUTINE ZPORFS1


      SUBROUTINE DPPTRS1( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: AP(*)
               REAL(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE DPPTRS
         END INTERFACE
         CALL DPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE DPPTRS1

      SUBROUTINE ZPPTRS1( UPLO, N, NRHS, AP, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: AP(*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZPPTRS
         END INTERFACE
         CALL ZPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE ZPPTRS1


      SUBROUTINE DPPRFS1( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, &
     &                    BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: AFP(*), AP(*), B(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX,  &
     &                         FERR, BERR, WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: AFP(*), AP(*), B( LDB,*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DPPRFS
         END INTERFACE
         REAL(WP) :: BERR1(1), FERR1(1)
         CALL DPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR1,    &
     &                BERR1, WORK, IWORK, INFO )
         BERR = BERR1(1); FERR = FERR1(1)
      END SUBROUTINE DPPRFS1

      SUBROUTINE ZPPRFS1( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR, &
     &                    BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX,  &
     &                         FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B( LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZPPRFS
         END INTERFACE
         REAL(WP) :: BERR1(1), FERR1(1)
         CALL ZPPRFS( UPLO, N, NRHS, AP, AFP, B, LDB, X, LDX, FERR1,    &
     &                BERR1, WORK, RWORK, INFO )
         BERR = BERR1(1); FERR = FERR1(1)
      END SUBROUTINE ZPPRFS1

      SUBROUTINE DPBTRS1( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AB( LDAB,*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB,     &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: AB( LDAB,*)
               REAL(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE DPBTRS
         END INTERFACE
         CALL DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
      END SUBROUTINE DPBTRS1

      SUBROUTINE ZPBTRS1( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB,     &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: AB( LDAB,*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZPBTRS
         END INTERFACE
         CALL ZPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
      END SUBROUTINE ZPBTRS1


      SUBROUTINE DPBRFS1( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B,   &
     &                    LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*), B(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, &
     &                         B, LDB, X, LDX, FERR, BERR, WORK, IWORK, &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N,    &
     &                                 NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*),     &
     &                                  B( LDB,*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DPBRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1) 
         CALL DPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB,  &
     &                X, LDX, FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DPBRFS1

      SUBROUTINE ZPBRFS1( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B,   &
     &                    LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, &
     &                         B, LDB, X, LDX, FERR, BERR, WORK, RWORK, &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) ::  KD, LDAB, LDAFB, LDB, LDX, N,    &
     &                                 NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) ::  AB( LDAB,*), AFB( LDAFB,*),  &
     &                                     B( LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZPBRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1) 
         CALL ZPBRFS( UPLO, N, KD, NRHS, AB, LDAB, AFB, LDAFB, B, LDB,  &
     &                X, LDX, FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZPBRFS1


      SUBROUTINE DPTTRS1( N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*)
         REAL(WP), INTENT(IN) :: E(*)
         REAL(WP), INTENT(OUT) :: B(*)
         INTERFACE
            SUBROUTINE DPTTRS( N, NRHS, D, E, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: D(*)
               REAL(WP), INTENT(IN) :: E(*)
               REAL(WP), INTENT(OUT) :: B( LDB,*)
            END SUBROUTINE DPTTRS
         END INTERFACE
         CALL DPTTRS( N, NRHS, D, E, B, LDB, INFO )
      END SUBROUTINE DPTTRS1

      SUBROUTINE ZPTTRS1( UPLO, N, NRHS, D, E, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*)
         COMPLEX(WP), INTENT(IN) :: E(*)
         COMPLEX(WP), INTENT(OUT) :: B(*)
         INTERFACE
            SUBROUTINE ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: D(*)
               COMPLEX(WP), INTENT(IN) :: E(*)
               COMPLEX(WP), INTENT(OUT) :: B( LDB,*)
            END SUBROUTINE ZPTTRS
         END INTERFACE
         CALL ZPTTRS( UPLO, N, NRHS, D, E, B, LDB, INFO )
      END SUBROUTINE ZPTTRS1


      SUBROUTINE DPTRFS1( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR,  &
     &                    BERR, WORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*), DF(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: B(*), E(*), EF(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX,   &
     &                         FERR, BERR, WORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: D(*), DF(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: B( LDB,*), E(*), EF(*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DPTRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL DPTRFS( N, NRHS, D, E, DF, EF, B, LDB, X, LDX, FERR1,     &
     &                BERR1, WORK, INFO )
      FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DPTRFS1

      SUBROUTINE ZPTRFS1( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,  &
     &                    FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: D(*), DF(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: B(*), E(*), EF(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X,  &
     &                         LDX, FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: D(*), DF(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: B( LDB,*), E(*), EF(*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZPTRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZPTRFS( UPLO, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,      &
     &                FERR1, BERR1, WORK, RWORK, INFO )
      FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZPTRFS1

      SUBROUTINE DSYTRS1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: A( LDA,*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,     &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER , INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(IN) :: A( LDA,*)
               REAL(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE DSYTRS
         END INTERFACE
         CALL DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      END SUBROUTINE DSYTRS1

      SUBROUTINE ZSYTRS1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,     &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER , INTENT(IN) :: IPIV(*)
               COMPLEX(WP), INTENT(IN) :: A( LDA,*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZSYTRS
         END INTERFACE
         CALL ZSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      END SUBROUTINE ZSYTRS1

      SUBROUTINE ZHETRS1( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER , INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: A( LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,     &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER , INTENT(IN) :: IPIV(*)
               COMPLEX(WP), INTENT(IN) :: A( LDA,*)
               COMPLEX(WP), INTENT(INOUT) :: B( LDB,*)
            END SUBROUTINE ZHETRS
         END INTERFACE
         CALL ZHETRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
      END SUBROUTINE ZHETRS1

      SUBROUTINE ZHERFS1( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,&
     &                    X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,&
     &                         LDB, X, LDX, FERR, BERR, WORK, RWORK,    &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*),      &
     &                                     B( LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZHERFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZHERFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, &
     &                LDX, FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZHERFS1

      SUBROUTINE ZSYRFS1( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,&
     &                    X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,&
     &                         LDB, X, LDX, FERR, BERR, WORK, RWORK,    &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*),      &
     &                                     B( LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X( LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZSYRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, &
     &                LDX, FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZSYRFS1

      SUBROUTINE DSYRFS1( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,&
     &                    X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*), B(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,&
     &                         LDB, X, LDX, FERR, BERR, WORK, IWORK,    &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDA, LDAF, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) ::  A( LDA,*), AF( LDAF,*),         &
     &                                  B( LDB,*)
               REAL(WP), INTENT(INOUT) :: X( LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DSYRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL DSYRFS( UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, &
     &                LDX, FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DSYRFS1

      SUBROUTINE DSPTRS1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(IN) :: AP(*)
               REAL(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE DSPTRS
         ENDINTERFACE
         CALL DSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE DSPTRS1

      SUBROUTINE ZSPTRS1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               COMPLEX(WP), INTENT(IN) :: AP(*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZSPTRS
         ENDINTERFACE
         CALL ZSPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE ZSPTRS1

      SUBROUTINE ZHPTRS1( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               COMPLEX(WP), INTENT(IN) :: AP(*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZHPTRS
         END INTERFACE
         CALL ZHPTRS( UPLO, N, NRHS, AP, IPIV, B, LDB, INFO )
      END SUBROUTINE ZHPTRS1

      SUBROUTINE ZHPRFS1( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, &
     &                    FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, &
     &                         LDX, FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZHPRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL ZHPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,     &
     &                FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZHPRFS1

      SUBROUTINE ZSPRFS1( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, &
     &                    FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, &
     &                         LDX, FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZSPRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL ZSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,     &
     &                FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZSPRFS1

      SUBROUTINE DSPRFS1( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX, &
     &                    FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         INTEGER, INTENT(IN) :: IPIV(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
         REAL(WP), INTENT(INOUT) :: X(LDX,*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, &
     &                         LDX, FERR, BERR, WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               INTEGER, INTENT(IN) :: IPIV(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: AFP(*), AP(*), B(LDB,*)
               REAL(WP), INTENT(INOUT) :: X(LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DSPRFS
         END INTERFACE
         REAL(WP) :: FERR1(1), BERR1(1)
         CALL DSPRFS( UPLO, N, NRHS, AP, AFP, IPIV, B, LDB, X, LDX,     &
     &                FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DSPRFS1

      SUBROUTINE DTRTRS1( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: A(LDA,*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B,   &
     &                         LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: A(LDA,*)
               REAL(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE DTRTRS
         END INTERFACE
         CALL DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,       &
     &                INFO )
      END SUBROUTINE DTRTRS1

      SUBROUTINE ZTRTRS1( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,   &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: A(LDA,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B,   &
     &                         LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: A(LDA,*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZTRTRS
         END INTERFACE
         CALL ZTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB,       &
     &                INFO )
      END SUBROUTINE ZTRTRS1

      SUBROUTINE DTRRFS1( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,&
     &                    LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: A(LDA,*), B(*)
         REAL(WP), INTENT(INOUT) :: X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B,   &
     &                         LDB, X, LDX, FERR, BERR, WORK, IWORK,    &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
               REAL(WP), INTENT(INOUT) :: X(LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DTRRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL DTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,    &
     &                LDX, FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DTRRFS1

      SUBROUTINE ZTRRFS1( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,&
     &                    LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(*)
         COMPLEX(WP), INTENT(INOUT) :: X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B,   &
     &                         LDB, X, LDX, FERR, BERR, WORK, RWORK,    &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDA, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: A(LDA,*), B(LDB,*)
               COMPLEX(WP), INTENT(INOUT) :: X(LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZTRRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZTRRFS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, X,    &
     &                LDX, FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZTRRFS1

      SUBROUTINE DTPTRS1( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,       &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AP(*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,  &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: AP(*)
               REAL(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE DTPTRS
         END INTERFACE
         CALL DTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE DTPTRS1

      SUBROUTINE ZTPTRS1( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,       &
     &                    INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AP(*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,  &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: AP(*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZTPTRS
         END INTERFACE
         CALL ZTPTRS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, INFO )
      END SUBROUTINE ZTPTRS1

      SUBROUTINE DTPRFS1( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X,    &
     &                    LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: AP(*), B(*), X(*)
         REAL(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,  &
     &                         X, LDX, FERR, BERR, WORK, IWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: AP(*), B(LDB,*), X(LDX,*)
               REAL(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE DTPRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL DTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,   &
     &                FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DTPRFS1

      SUBROUTINE ZTPRFS1( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X,    &
     &                    LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AP(*), B(*), X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB,  &
     &                         X, LDX, FERR, BERR, WORK, RWORK, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AP(*), B(LDB,*), X(LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZTPRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZTPRFS( UPLO, TRANS, DIAG, N, NRHS, AP, B, LDB, X, LDX,   &
     &                FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZTPRFS1

      SUBROUTINE DTBTRS1( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,  &
     &                    LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(IN) :: AB(LDAB,*)
         REAL(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB,&
     &                         B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(IN) :: AB(LDAB,*)
               REAL(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE DTBTRS
         END INTERFACE
         CALL DTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, &
     &                INFO )
      END SUBROUTINE DTBTRS1

      SUBROUTINE ZTBTRS1( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,  &
     &                    LDB, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         COMPLEX(WP), INTENT(IN) :: AB(LDAB,*)
         COMPLEX(WP), INTENT(INOUT) :: B(*)
         INTERFACE
            SUBROUTINE ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB,&
     &                         B, LDB, INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               COMPLEX(WP), INTENT(IN) :: AB(LDAB,*)
               COMPLEX(WP), INTENT(INOUT) :: B(LDB,*)
            END SUBROUTINE ZTBTRS
         END INTERFACE
         CALL ZTBTRS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, &
     &                INFO )
      END SUBROUTINE ZTBTRS1

      SUBROUTINE DTBRFS1( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,  &
     &                    LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO, IWORK(*)
         REAL(WP), INTENT(OUT) :: BERR, FERR
         REAL(WP), INTENT(IN) :: AB(LDAB,*), B(*), X(*)
         REAL(WP), INTENT(IN) :: WORK(*)
         INTERFACE
            SUBROUTINE DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB,&
     &                         B, LDB, X, LDX, FERR, BERR, WORK, IWORK, &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO, IWORK(*)
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*)
               REAL(WP), INTENT(IN) :: AB(LDAB,*), B(LDB,*), X(LDX,*)
               REAL(WP), INTENT(IN) :: WORK(*)
            END SUBROUTINE DTBRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL DTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, &
     &                X, LDX, FERR1, BERR1, WORK, IWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE DTBRFS1

      SUBROUTINE ZTBRFS1( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B,  &
     &                    LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO )
         USE LA_PRECISION, ONLY: WP => DP
         CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
         INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
         INTEGER, INTENT(OUT) :: INFO
         REAL(WP), INTENT(OUT) :: BERR, FERR, RWORK(*)
         COMPLEX(WP), INTENT(IN) :: AB(LDAB,*), B(*), X(*)
         COMPLEX(WP), INTENT(OUT) :: WORK(*)
         INTERFACE
            SUBROUTINE ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB,&
     &                         B, LDB, X, LDX, FERR, BERR, WORK, RWORK, &
     &                         INFO )
               USE LA_PRECISION, ONLY: WP => DP
               CHARACTER(LEN=1), INTENT(IN) :: DIAG, TRANS, UPLO
               INTEGER, INTENT(IN) :: KD, LDAB, LDB, LDX, N, NRHS
               INTEGER, INTENT(OUT) :: INFO
               REAL(WP), INTENT(OUT) :: BERR(*), FERR(*), RWORK(*)
               COMPLEX(WP), INTENT(IN) :: AB(LDAB,*), B(LDB,*), X(LDX,*)
               COMPLEX(WP), INTENT(OUT) :: WORK(*)
            END SUBROUTINE ZTBRFS
         END INTERFACE
         REAL(WP) FERR1(1), BERR1(1)
         CALL ZTBRFS( UPLO, TRANS, DIAG, N, KD, NRHS, AB, LDAB, B, LDB, &
     &                X, LDX, FERR1, BERR1, WORK, RWORK, INFO )
         FERR = FERR1(1); BERR = BERR1(1)
      END SUBROUTINE ZTBRFS1



      END MODULE F77_LAPACK
