!****h* simple/invert_gpu_utils
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME
! invert_gpu_utils - Various utilities for matrix inversion double
!
! AUTHOR
! Frederic D.R. Bonnet.
!
! DESCRIPTION
! Various utilities for the simple_math_gpu modules:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
module invert_cpu_utils

  use simple_defs
  use simple_timing

  implicit none
 
  logical,external :: LSAME

!******

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the LU decomposition of a Double precision 2d array of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getLU_SMatrix_cpu(n,m,lda,size_of_elt,lwork,matrix)
    implicit none
    integer                         :: n,m
    integer                         :: lda
    integer                         :: size_of_elt
    integer                         :: lwork

    real(sp)                        :: matrix(LDA,*)

    !local variable

    integer                         :: err
    integer,allocatable             :: ipiv(:)
    real(sp),allocatable            :: work(:)

    allocate(ipiv(n))
    allocate(work(lwork))

    call sgetrf(n, m, matrix, n, ipiv, err )!my wrapper call to magma

    deallocate(ipiv)
    deallocate(work)

    return
  end subroutine getLU_SMatrix_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the LU decomposition of a Double precision 2d array of size
! (n)x(m) on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getLU_DMatrix_cpu(n,m,lda,size_of_elt,lwork,matrix)
    implicit none
    integer                         :: n,m
    integer                         :: lda
    integer                         :: size_of_elt
    integer                         :: lwork

    real(dp)                        :: matrix(LDA,*)

    !local variable

    integer                         :: err
    integer,allocatable             :: ipiv(:)
    real(dp),allocatable            :: work(:)

    allocate(ipiv(n))
    allocate(work(lwork))

    call dgetrf(n, m, matrix, n, ipiv, err )!my wrapper call to magma

    deallocate(ipiv)
    deallocate(work)

    return
  end subroutine getLU_DMatrix_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a Double precision 2d array of size (n)x(m)
! on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseSMatrix_cpu(n,m,lda,size_of_elt,lwork,A, invSA)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    real(sp)                      :: invSA(LDA,*)
    real(sp)                      :: A(LDA,*)

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    real(sp),allocatable          :: work(:)

    !timer variables

    double precision              :: elps_sgetrf
    double precision              :: elps_sgetri
    double precision,dimension(2) :: st_sgetrf
    double precision,dimension(2) :: et_sgetrf
    double precision,dimension(2) :: st_sgetri
    double precision,dimension(2) :: et_sgetri

    !start of the execution commands

    allocate(ipiv(n))
    allocate(work(lwork))

    invSA(1:n,1:m) = A(1:n,1:m)
    
#if defined (BENCH) && defined (LINUX)
    call gettimeofday_c(st_sgetrf)
    call start_timer_cpu("sgetrf")
#endif
    call sgetrf(n, m, invSA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "the LU decomposition in sgetrf has failed"
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("sgetrf")
    call gettimeofday_c(et_sgetrf)
    call elapsed_time_c(st_sgetrf,et_sgetrf,elps_sgetrf)

    call gettimeofday_c(st_sgetri)
    call start_timer_cpu("sgetri")
#endif
    call sgetri( n, invSA, lda, ipiv, work, lwork, err )
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("sgetri")
    call gettimeofday_c(et_sgetri)
    call elapsed_time_c(st_sgetri,et_sgetri,elps_sgetri)
#endif

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH) && defined (LINUX)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_sgetrf, ", ", &
         elps_sgetri, ", ", &
         (elps_sgetrf + elps_sgetri), ", ", &
         Log10(elps_sgetrf), ", ", &
         Log10(elps_sgetri), ", ", &
         Log10(elps_sgetrf+elps_sgetri)
    !end of benchmark output
#endif

    return
  end subroutine getBlockInverseSMatrix_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a Double precision 2d array of size (n)x(m)
! on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseDMatrix_cpu(n,m,lda,size_of_elt,lwork,A, invDA)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    real(dp)                      :: invDA(LDA,*)
    real(dp)                      :: A(LDA,*)

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    real(dp),allocatable          :: work(:)

    !timer variables

    double precision              :: elps_dgetrf
    double precision              :: elps_dgetri
    double precision,dimension(2) :: st_dgetrf
    double precision,dimension(2) :: et_dgetrf
    double precision,dimension(2) :: st_dgetri
    double precision,dimension(2) :: et_dgetri

    !start of the execution commands

    allocate(ipiv(n))
    allocate(work(lwork))

    invDA(1:n,1:m) = A(1:n,1:m)
    
#if defined (BENCH) && defined (LINUX)
    call gettimeofday_c(st_dgetrf)
    call start_timer_cpu("dgetrf")
#endif
    call dgetrf(n, m, invDA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "the LU decomposition in dgetrf has failed"
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("dgetrf")
    call gettimeofday_c(et_dgetrf)
    call elapsed_time_c(st_dgetrf,et_dgetrf,elps_dgetrf)

    call gettimeofday_c(st_dgetri)
    call start_timer_cpu("dgetri")
#endif
    call dgetri( n, invDA, lda, ipiv, work, lwork, err )
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("dgetri")
    call gettimeofday_c(et_dgetri)
    call elapsed_time_c(st_dgetri,et_dgetri,elps_dgetri)
#endif

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH) && defined (LINUX)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_dgetrf, ", ", &
         elps_dgetri, ", ", &
         (elps_dgetrf + elps_dgetri), ", ", &
         Log10(elps_dgetrf), ", ", &
         Log10(elps_dgetri), ", ", &
         Log10(elps_dgetrf+elps_dgetri)
    !end of benchmark output
#endif

    return
  end subroutine getBlockInverseDMatrix_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a Double complex precision 2d array of
! size (n)x(m)
! on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseCMatrix_cpu(n,m,lda,size_of_elt,lwork,A, invCA)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    complex(sp)                   :: invCA(LDA,*)
    complex(sp)                   :: A(LDA,*)

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    complex(sp),allocatable       :: work(:)

    !timer variables

    double precision              :: elps_cgetrf
    double precision              :: elps_cgetri
    double precision,dimension(2) :: st_cgetrf
    double precision,dimension(2) :: et_cgetrf
    double precision,dimension(2) :: st_cgetri
    double precision,dimension(2) :: et_cgetri

    !start of the execution commands

    allocate(ipiv(n))
    allocate(work(lwork))

    invCA(1:n,1:m) = A(1:n,1:m)
    
#if defined (BENCH) && defined (LINUX)
    call gettimeofday_c(st_cgetrf)
    call start_timer_cpu("cgetrf")
#endif
    call cgetrf(n, m, invCA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "the LU decomposition in cgetrf has failed"
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("cgetrf")
    call gettimeofday_c(et_cgetrf)
    call elapsed_time_c(st_cgetrf,et_cgetrf,elps_cgetrf)

    call gettimeofday_c(st_cgetri)
    call start_timer_cpu("cgetri")
#endif
    call cgetri( n, invCA, lda, ipiv, work, lwork, err )
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("cgetri")
    call gettimeofday_c(et_cgetri)
    call elapsed_time_c(st_cgetri,et_cgetri,elps_cgetri)
#endif

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH) && defined (LINUX)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_cgetrf, ", ", &
         elps_cgetri, ", ", &
         (elps_cgetrf + elps_cgetri), ", ", &
         Log10(elps_cgetrf), ", ", &
         Log10(elps_cgetri), ", ", &
         Log10(elps_cgetrf+elps_cgetri)
    !end of benchmark output
#endif

    return
  end subroutine getBlockInverseCMatrix_cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a Double complex precision 2d array of
! size (n)x(m)
! on CPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseZMatrix_cpu(n,m,lda,size_of_elt,lwork,A, invZA)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    complex(dp)                   :: invZA(LDA,*)
    complex(dp)                   :: A(LDA,*)

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    complex(dp),allocatable       :: work(:)

    !timer variables

    double precision              :: elps_zgetrf
    double precision              :: elps_zgetri
    double precision,dimension(2) :: st_zgetrf
    double precision,dimension(2) :: et_zgetrf
    double precision,dimension(2) :: st_zgetri
    double precision,dimension(2) :: et_zgetri

    !start of the execution commands

    allocate(ipiv(n))
    allocate(work(lwork))

    invZA(1:n,1:m) = A(1:n,1:m)
    
#if defined (BENCH) && defined (LINUX)
    call gettimeofday_c(st_zgetrf)
    call start_timer_cpu("zgetrf")
#endif
    call zgetrf(n, m, invZA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "the LU decomposition in zgetrf has failed"
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("zgetrf")
    call gettimeofday_c(et_zgetrf)
    call elapsed_time_c(st_zgetrf,et_zgetrf,elps_zgetrf)

    call gettimeofday_c(st_zgetri)
    call start_timer_cpu("zgetri")
#endif
    call zgetri( n, invZA, lda, ipiv, work, lwork, err )
#if defined (BENCH) && defined (LINUX)
    call stop_timer_cpu("zgetri")
    call gettimeofday_c(et_zgetri)
    call elapsed_time_c(st_zgetri,et_zgetri,elps_zgetri)
#endif

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH) && defined (LINUX)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_zgetrf, ", ", &
         elps_zgetri, ", ", &
         (elps_zgetrf + elps_zgetri), ", ", &
         Log10(elps_zgetrf), ", ", &
         Log10(elps_zgetri), ", ", &
         Log10(elps_zgetrf+elps_zgetri)
    !end of benchmark output
#endif

    return
  end subroutine getBlockInverseZMatrix_cpu

end module invert_cpu_utils
