! ============================================================================
! Name        : testing_SMatvect_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 24th of March 2016
! Description : tests the sgemv within the CPU environment
!             :
! ============================================================================
!
program testing_SMatvect_cpu
  use simple_defs
  use simple_cuda_defs
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu
  use simple_math, only: calc_corr, csq
  use simple_yaml_output
  use simple_yaml_strings
  !packages needed if in the simple{_gpu,_test} environment to include in *.csh
!  use invert_cpu_utils  
!  use matmul_cpu_utils
!  use simple_cuda_defs
!  use matmul_gpu_utils
!  use invert_gpu_utils
  !packages needed needed if the testcode environment other already in *.csh
  use simple_math_gpu

  !$ use omp_lib
  !$ use omp_lib_kinds

  implicit none
  logical, parameter :: ibench_local=.true.        !< benchmark result or not
  logical, parameter :: ibench_write_local=.true.  !< write benchmark result or not

  logical, parameter :: debug=.true.          !< debug indicator
  logical, parameter :: debug_cpu=.false.     !false/true
  logical, parameter :: debug_high=.true.     !true/false
  logical, parameter :: debug_write=.false.   !false/true
  logical, parameter :: debug_write_C=.false. !false/true
  
  logical, parameter :: call_cpu_corr=.true.!< debug indicator
  
#define devptr_t integer*8

  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
  integer, parameter            :: npart= 35    !4000  !30000
  integer, parameter            :: start_npart=35, istep=64
  integer, parameter            :: nrot = 10000 !4000 !15000
  
  integer                       :: omp_nthr=8           !openMP variables

  !gpu gear
  integer                       :: ndev
  integer                       :: rc
  integer                       :: size_pft1, size_pft2
  integer                       :: size_prod
  integer                       :: err
  integer                       :: lda,incx,incy
  real(sp)                      :: alpha,beta
  real(dp)                      :: gflops

  !comparators
  real(sp)                      :: sum_diff_prod
  !variables for the function 
  real(sp), allocatable         :: pft1(:,:)
  real(sp), allocatable         :: pft2(:)
  real(sp), allocatable         :: prod_cpu(:)
  real(sp), allocatable         :: prod_gpu(:)
  !timer variables

  real(dp)                      :: speedup

  double precision              :: elps_mtml_gpu
  double precision,dimension(2) :: st_ZM_r, et_ZM_r

  double precision              :: elps_mtml_cpu
  double precision,dimension(2) :: st_mtml_cpu, et_mtml_cpu

  !indexers
  integer                       :: i,ik,irot
  integer                       :: idev
  !index variable major
  integer                       :: ipart

  !functions calls
  integer :: get_dev_count_c
  
  !start of the execution commandsinterface

  !Start of the OpenMP environment
  !$ call omp_set_num_threads(omp_nthr)

  !first setup the debugging options first using 0 and 1.
  if (        debug) s_debug_gpu%debug_i         = 1 !else 0
  if (    debug_cpu) s_debug_gpu%debug_cpu_i     = 1
  if (   debug_high) s_debug_gpu%debug_high_i    = 1
  if (  debug_write) s_debug_gpu%debug_write_i   = 1
  if (debug_write_C) s_debug_gpu%debug_write_C_i = 1
#if defined (BENCH)
  ibench = .true.
  ibench_write = .true.
  if (ibench) s_bench%bench_i = 1
  if (ibench_write) s_bench%bench_write_i = 1
#endif
  !start of the greeting message
  call timestamp()
  call start_Alltimers_cpu()

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  call Sanity_check_cpu(sysQ, hstD)

  alpha = 1.0
  beta = 0.0
  do ipart = start_npart, npart, istep

     lda = ipart
     incx = 1
     incy = 1
     
     gflops = dble(ipart*nrot*1.0d0)/1.e9

!*******************************************************************************
!    Data allocation and initialisation of the matrices for testing purposes.
!
!*******************************************************************************

     !cpu allocation
     allocate(pft1(ipart, nrot)) !pft1
     allocate(pft2(       nrot)) !pft2

     allocate(prod_cpu(nrot)) !product matrix
     allocate(prod_gpu(nrot)) !product matrix

     pft1 = 0.0
     pft2 = 0.0
     prod_cpu = 0.0
     prod_gpu = 0.0
     
!*******************************************************************************
!    Filling and printing PFTs array with random numbers using test functions
!
!*******************************************************************************
     call getSRandomGaussianDistr_2D(ipart, nrot, pft1)
     call getSRandomGaussianDistr_1D(       nrot, pft2)

     write(*,*)
     write(*,'(10x,a,5x,a,6x,a)')"ipart","irot","(pft1)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
              write(*,*)i, irot,pft1(i,irot)
        end do
     end do
     write(*,*)
     write(*,'(5x,a,6x,a)')"irot","(pft2)"
     do irot=1,nrot - (nrot-4)
        write(*,*) irot, pft2(irot)
     end do

!*******************************************************************************
!     now testign the mtml function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************
     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     Now testing the matmul versus gemvfunction                '
     write(*,*)'***************************************************************'
     size_pft1 = ipart*nrot
     size_pft2 = nrot
     size_prod = nrot
     write(*,*)
     write(*,'(3x,a,3x,a,4x,a,2x,a,2x,a)') &
          "npart","ipart","nrot","(1*ipart)","OpenMP thrd"
     write(*,'(8(3x,i5))') npart, ipart, nrot, (1*ipart), omp_nthr
     write(*,*) "N (Ref)       pft1(npart,nrot): ",size_pft1
     write(*,*) "N (Ref)       pft1(npart,nrot): ",size_pft1/1.e6,"(Millions)"
     write(*,*) "N (ptc)       pft2(nrot)      : ",size_pft2
     write(*,*) "N (ptc)       pft2(nrot)      : ",size_pft2/1.e6,"(Millions)"
     write(*,*) " Ratio of size (PFT2 and PFT1): ",size_pft2 / real(size_pft1)
     write(*,*)'                                                               '
     write(*,*)'***************CPU Matrix Mul with matvec (OpenMP)*************'
     write(*,*)'                                                               '
     
     call start_timer_cpu("matvec_cpu")
     call gettimeofday_c(st_mtml_cpu)
     
!     !$omp parallel default(shared)
!     !$omp workshare
     prod_cpu = matmul(pft1,pft2)
!     !$omp end workshare nowait
!     !$omp end parallel
     
     call gettimeofday_c(et_mtml_cpu)
     call elapsed_time_c(st_mtml_cpu,et_mtml_cpu,elps_mtml_cpu)
     call stop_timer_cpu("matvec_cpu")

     write(*,*)"Elapsed time matvec(OMP)  : ",real(elps_mtml_cpu),"(seconds)"
     write(*,*)'                                                               '
     write(*,*)'***************CPU Matrix Mul with sgemv***********************'
     write(*,*)'                                                               '

     call start_timer_cpu("mtml_gpu")
     call gettimeofday_c(st_ZM_r)

     call sgemv("N", ipart, nrot, &
                alpha,            &
                pft1, lda,        &
                pft2, incx,       &
                beta,             &
                prod_gpu,incy     )
     call gettimeofday_c(et_ZM_r)
     call elapsed_time_c(st_ZM_r,et_ZM_r,elps_mtml_gpu)
     call stop_timer_cpu("mtml_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_polarft_multi_gpus_gpu_c_ = -1"
     
     write(*,*)"Elapsed time Sgemv        : ",real(elps_mtml_gpu),"(seconds)"
     write(*,*)'                                                               '
     write(*,*)'*********Comparing product from matvec vs Sgemv (CPU)**********'
     write(*,*)'                                                               '

     write(*,'(5x,a,5x,a,4x,a)') &
          "irot","(prod_cpu)","(prod_gpu)"
     do irot=1,nrot - (nrot-5)
        write(*,*)irot, prod_cpu(irot), prod_gpu(irot)
     end do

     sum_diff_prod = sum(prod_cpu-prod_gpu)
     write(*,*) "Diff sum(prod_matmul-prod_sgemv): ",sum_diff_prod
     write(*,*) "N           prod_{cpu,gpu}(nrot): ",size_prod/1.e6,"(Millions)!"
     write(*,*)'***************************************************************'
     speedup = elps_mtml_cpu / elps_mtml_gpu
     gflops = gflops / elps_mtml_gpu
     write(*,'(3x,a,2x,a,6x,a)') "Matvec CPU (secs)","Sgemm CPU (secs)","Gflp/s"
     write(*,'(x,f15.5,3x,f15.5,7x,f10.5)')elps_mtml_cpu,elps_mtml_gpu,gflops
     write(*,*)
     write(*,'(x,a,x,f10.3)')" Speed up from matmul vs sgemv : ",speedup
     write(*,*)'***************************************************************'

#if defined (BENCH)
     if( ibench_write ) then
        open(4,file='SMatvect_2D.asc',status='unknown',position='append')
        write(4,'(1x,f15.8,2x,f15.8)')ipart*nrot/1.e6, &
             real(elps_mtml_cpu) / real(elps_mtml_gpu)
        close(4)
     end if
#endif     

!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************

     !Freeing the CPU memory 
     deallocate(pft1)
     deallocate(pft2)
     deallocate(prod_cpu)
     deallocate(prod_gpu)

  end do !end of do ipart = start_npart, npart, istep
  
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
  

  !shutting down the timers
  call stop_Alltimers_cpu()

contains

end program testing_SMatvect_cpu
!*******************************************************************************
!    Subroutine to run sanity checks on the data structure passed CPU
!
!*******************************************************************************
!
subroutine Sanity_check_cpu(sysQ, hstD)
  use simple_defs
  use greeting_version
  use simple_systemQuery_cpu
  implicit none
  type(systemQuery_cpu)           :: sysQ
  type(systemDetails)             :: hstD
  !local variables
  !cpu gear
  integer                         :: nCPU_cores
  integer*8                       :: h_TotalMemSize
  integer*8                       :: h_AvailMemSize

  !start of the execution commands

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'  Sanity checks on the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  h_TotalMemSize = sysQ%get_SysTotalMem_size()
  h_AvailMemSize = sysQ%get_SysAvailMem_size()
  write(*,*)'Number of cores on system     (sysQ):',nCPU_cores
  write(*,*)'Number of cores               (hstD):',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'Total Mem on system           (sysQ):',h_TotalMemSize
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_Size
  if(h_TotalMemSize /= hstD%mem_Size)call sysQ%get_warning_dataStruct()
  write(*,*)'Total Available Mem on system (sysQ):',h_AvailMemSize
#if defined (MACOSX)
  write(*,*)'Total Mem on system           (hstD):',hstD%mem_User
  if(h_AvailMemSize /= hstD%mem_User)call sysQ%get_warning_dataStruct()
#elif defined (LINUX)
  write(*,*)'Total Mem on system           (hstD):',hstD%avail_Mem
  if(h_AvailMemSize /= hstD%avail_Mem)call sysQ%get_warning_dataStruct()
#endif
  
  return
end subroutine Sanity_check_cpu
