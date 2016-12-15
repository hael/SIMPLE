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
module invert_gpu_utils

  use simple_defs
  use simple_cuda_defs
  use simple_timing

  implicit none

#define size_t integer*8
#define devptr_t integer*8
  
  logical,external :: LSAME

!******

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the LU decomposition of a Double precision 2d array of size
! (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getLU_DMatrix_gpu(n,m,lda,size_of_elt,lwork,devPtrA,LU_matrix)
    implicit none
    integer                         :: n,m
    integer                         :: lda
    integer                         :: size_of_elt
    integer                         :: lwork

    real(dp)                        :: LU_matrix(LDA,*)

    devptr_t                        :: devPtrA

    !local variable

    integer                         :: err
    integer,allocatable             :: ipiv(:)
    real(dp),allocatable            :: work(:)

#if defined (CUDA) && defined (MAGMA)

    allocate(ipiv(n))
    allocate(work(lwork))

    call magmaf_init()
    call magma_dgetrf_gpu(n, m, devPtrA, n, ipiv, err )!my wrapper call to magma
    err = cublas_get_matrix(n, m, size_of_elt, devPtrA, m, LU_matrix, m)

    deallocate(ipiv)
    deallocate(work)

    !finalising the magma environment
    call magmaf_finalize()

#endif

    return
  end subroutine getLU_DMatrix_gpu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a Double precision 2d array of size (n)x(m)
! on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseDMatrix_gpu(n,m,lda,size_of_elt,lwork,devPtrA, &
                                        invDA_gpu)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    real(dp)                      :: invDA_gpu(LDA,*)

    devptr_t                      :: devPtrA

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    real(dp),allocatable          :: work(:)

    !timer variables

    double precision              :: elps_mgblk_dgetrf
    double precision              :: elps_mgblk_dgetri
    double precision,dimension(2) :: st_mgblk_dgetrf
    double precision,dimension(2) :: et_mgblk_dgetrf
    double precision,dimension(2) :: st_mgblk_dgetri
    double precision,dimension(2) :: et_mgblk_dgetri

    !start of the execution commands
    
#if defined (MAGMA)

!    write(*,*)" hello subroutine getBlockInverseDMatrix_gpu"

    allocate(ipiv(n))
    allocate(work(lwork))

    call magmaf_init()

#if defined (BENCH)
    call gettimeofday_c(st_mgblk_dgetrf)
    call start_timer_cpu("magma_dgetrf_gpu")
#endif
    call magma_dgetrf_gpu(n, m, devPtrA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "the LU decomposition in magma_dgetrf_gpu has failed"
#if defined (BENCH)
    call stop_timer_cpu("magma_dgetrf_gpu")
    call gettimeofday_c(et_mgblk_dgetrf)
    call elapsed_time_c(st_mgblk_dgetrf,et_mgblk_dgetrf,elps_mgblk_dgetrf)

    call gettimeofday_c(st_mgblk_dgetri)
    call start_timer_cpu("magma_dgetri_gpu")
#endif
    call magma_dgetri_block_gpu_v2( n, devPtrA, lda, ipiv, work, lwork, err )
#if defined (BENCH)
    call stop_timer_cpu("magma_dgetri_gpu")
    call gettimeofday_c(et_mgblk_dgetri)
    call elapsed_time_c(st_mgblk_dgetri,et_mgblk_dgetri,elps_mgblk_dgetri)
#endif

    err = cublas_get_matrix( n, m, size_of_elt, devPtrA, m, invDA_gpu, m)

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_mgblk_dgetrf, ", ", &
         elps_mgblk_dgetri, ", ", &
         (elps_mgblk_dgetrf + elps_mgblk_dgetri), ", ", &
         Log10(elps_mgblk_dgetrf), ", ", &
         Log10(elps_mgblk_dgetri), ", ", &
         Log10(elps_mgblk_dgetrf+elps_mgblk_dgetri)
    !end of benchmark output
#endif
    !finalising the magma environment
    call magmaf_finalize()
#endif

    return
  end subroutine getBlockInverseDMatrix_gpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION
! subroutine to take the inverse of a DoubleComplex precision 2d array of
! size (n)x(m) on GPU 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SOURCE
  subroutine getBlockInverseZMatrix_gpu(n,m,lda,size_of_elt,lwork,devPtrA,&
                                        invZA_gpu)
    implicit none

    !global variables

    integer                       :: n,m
    integer                       :: lda
    integer                       :: size_of_elt
    integer                       :: lwork

    complex(dp)                   :: invZA_gpu(LDA,*)

    devptr_t                      :: devPtrA

    !local variable

    integer                       :: err
    integer,allocatable           :: ipiv(:)
    complex(dp),allocatable       :: work(:)

    !timer variables

    double precision              :: elps_mgblk_zgetrf
    double precision              :: elps_mgblk_zgetri
    double precision,dimension(2) :: st_mgblk_zgetrf
    double precision,dimension(2) :: et_mgblk_zgetrf
    double precision,dimension(2) :: st_mgblk_zgetri
    double precision,dimension(2) :: et_mgblk_zgetri

    !start of the execution commands
    
#if defined (MAGMA)

    allocate(ipiv(n))
    allocate(work(lwork))

    !initialising the magma environment
    call magmaf_init()

#if defined (BENCH)
    call gettimeofday_c(st_mgblk_zgetrf)
    call start_timer_cpu("magma_zgetrf_gpu")
#endif
    call magma_zgetrf_gpu(n, m, devPtrA, n, ipiv, err )!my wrapper call to magma
    if (err .ne. 0) write(*,*) "LU decomposition in magma_zgetrf_gpu has failed"
#if defined (BENCH)
    call stop_timer_cpu("magma_zgetrf_gpu")
    call gettimeofday_c(et_mgblk_zgetrf)
    call elapsed_time_c(st_mgblk_zgetrf,et_mgblk_zgetrf,elps_mgblk_zgetrf)

    call gettimeofday_c(st_mgblk_zgetri)
    call start_timer_cpu("magma_zgetri_gpu")
#endif
    call magma_zgetri_block_gpu_v2( n, devPtrA, lda, ipiv, work, lwork, err )
#if defined (BENCH)
    call stop_timer_cpu("magma_zgetri_gpu")
    call gettimeofday_c(et_mgblk_zgetri)
    call elapsed_time_c(st_mgblk_zgetri,et_mgblk_zgetri,elps_mgblk_zgetri)
#endif

    err = cublas_get_matrix(n, m, size_of_elt, devPtrA, m, invZA_gpu, m)

    deallocate(ipiv)
    deallocate(work)

#if defined (BENCH)
    !benchmark output.
    write(2500,'(i5,a,f10.5,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8,a,f15.8)') &
         n, ", ", &
         log10(dble(n)), ", ", &
         elps_mgblk_zgetrf, ", ", &
         elps_mgblk_zgetri, ", ", &
         (elps_mgblk_zgetrf + elps_mgblk_zgetri), ", ", &
         Log10(elps_mgblk_zgetrf), ", ", &
         Log10(elps_mgblk_zgetri), ", ", &
         Log10(elps_mgblk_zgetrf+elps_mgblk_zgetri)
    !end of benchmark output
#endif

    !finalising the magma environment
    call magmaf_finalize()

!    write(*,*)" bye subroutine getBlockInverseZMatrix_gpu"

#endif

    return
  end subroutine getBlockInverseZMatrix_gpu

end module invert_gpu_utils
