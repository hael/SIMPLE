! ============================================================================
! Name        : testing_SMatmult_cpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 21th of March 2016
! Description : tests the sgemm within the CPU environment
!             :
! ============================================================================
!
program testing_SMatmult_gpu
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_deviceQuery_gpu
  use simple_systemQuery_cpu
  use simple_math, only: calc_corr, csq
  use simple_yaml_output
  use simple_yaml_strings
  !packages needed if in the simple{_gpu,_test} environment to include in *.csh
  !use invert_cpu_utils  
  !use matmul_cpu_utils
  !use simple_cuda_defs
  !use matmul_gpu_utils
  !use invert_gpu_utils
  !packages needed needed if the testcode environment other already in *.csh
  use simple_math_gpu

  !$ use omp_lib
  !$ use omp_lib_kinds

  implicit none
  logical, parameter :: ibench_local=.true.       !benchmark result or not
  logical, parameter :: ibench_write_local=.true. !write benchmark result or not

  logical, parameter :: debug=.true.          !< debug indicator
  logical, parameter :: debug_cpu=.false.     !false/true
  logical, parameter :: debug_high=.true.     !true/false
  logical, parameter :: debug_write=.false.   !false/true
  logical, parameter :: debug_write_C=.false. !false/true
  
  logical, parameter :: call_cpu_corr=.true.!< debug indicator
  
#define devptr_t integer*8

  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
  integer, parameter            :: npart= 4000   !4000  !30000
  integer, parameter            :: start_npart=4000, istep=64
  integer, parameter            :: nrot = 4000 !4000 !15000
  integer, parameter            :: nk = 1000
  
  integer                       :: omp_nthr=8           !openMP variables

  !local variables
  type(deviceDetails),allocatable :: a_devD(:)
  !gpu gear
  integer                       :: ndev
  integer                       :: rc
  integer                       :: size_pft1, size_pft2
  integer                       :: size_prod
  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc
  real(sp)                      :: alpha,beta
  real(dp)                      :: gflops

  !comparators
  real(sp)                      :: sum_diff_prod
  !variables for the function 
  real(sp), allocatable         :: pft1(:,:)
  real(sp), allocatable         :: pft2(:,:)
  real(sp), allocatable         :: prod_cpu(:,:)
  real(sp), allocatable         :: prod_gpu(:,:)
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
  call hello_deviceQuery_gpu(err)
  call timestamp()
  call start_Alltimers_cpu()

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  call Sanity_check_cpu(sysQ, hstD)
  write(*,*)'******************************************************************'
  write(*,*)'   Device fills in the object(devQ) and the data structure(devD)  '
  write(*,*)'******************************************************************'
  rc = get_dev_count_c(devD)
  ndev = devD%ndev
  allocate(a_devD(0:devD%ndev-1))
  do idev = 0, ndev-1
     call devQ%new_deviceQuery_gpu(devD,idev)
     !call Sanity_check_gpu(devQ, devD)
     !mapping the data structures into an array
     a_devD(idev) = devD
  end do

  call Sanity_check_a_devD_gpu(a_devD,ndev)

  alpha = 1.0
  beta = 0.0
  do ipart = start_npart, npart, istep

     lda = ipart
     ldb = nk
     ldc = ipart

     gflops = dble(ipart*nrot*1.0d0)/1.e9
     gflops = gflops * dble(nk)

!*******************************************************************************
!    Data allocation and initialisation of the matrices for testing purposes.
!
!*******************************************************************************

     !cpu allocation
     allocate(pft1(ipart, nk)) !pft1
     allocate(pft2( nk, nrot)) !pft2

     allocate(prod_cpu(ipart, nrot)) !product matrix
     allocate(prod_gpu(ipart, nrot)) !product matrix

     pft1 = 0.0
     pft2 = 0.0
     prod_cpu = 0.0
     prod_gpu = 0.0
     
!*******************************************************************************
!    Filling and printing PFTs array with random numbers using test functions
!
!*******************************************************************************
     call getSRandomGaussianDistr_2D(ipart,   nk, pft1)
     call getSRandomGaussianDistr_2D(   nk, nrot, pft2)

     write(*,*)
     write(*,'(10x,a,5x,a,6x,a,10x,a)')"ipart","ik","(pft1)","(pft2)"
     do i=1,ipart - (ipart-2)
        do ik=1,nk - (nk-2)
              write(*,*)i, ik,pft1(i,ik), pft2(ik,i)
        end do
     end do
!*******************************************************************************
!     now testign the mtml function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************
     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     Now testing the mtml function                             '
     write(*,*)'***************************************************************'
     size_pft1 = ipart*nrot
     size_pft2 = ipart*nrot
     size_prod = ipart*ipart
     write(*,*)
     write(*,'(3x,a,3x,a,4x,a,2x,a,2x,a)') &
          "npart","ipart","nrot","(1*ipart)","OpenMP thrd"
     write(*,'(8(3x,i5))') npart, ipart, nrot, (1*ipart), omp_nthr
     write(*,*) "N (Ref)       pft1(npart,nrot): ",size_pft1
     write(*,*) "N (Ref)       pft1(npart,nrot): ",size_pft1/1.e6,"(Millions)"
     write(*,*) "N (ptc)       pft2(npart,nrot): ",size_pft2
     write(*,*) "N (ptc)       pft2(npart,nrot): ",size_pft2/1.e6,"(Millions)"
     write(*,*) " Ratio of size (PFT2 and PFT1): ",size_pft2 / real(size_pft1)
     write(*,*)'                                                               '
     write(*,*)'***************CPU Matrix Mul with matmul (OpenMP)*************'
     write(*,*)'                                                               '
     
     call start_timer_cpu("matmul_cpu")
     call gettimeofday_c(st_mtml_cpu)
     
     !$omp parallel default(shared)
     !$omp workshare
     prod_cpu = matmul(pft1,pft2)
     !$omp end workshare nowait
     !$omp end parallel
     
     call gettimeofday_c(et_mtml_cpu)
     call elapsed_time_c(st_mtml_cpu,et_mtml_cpu,elps_mtml_cpu)
     call stop_timer_cpu("matmul_cpu")

     write(*,*)"Elapsed time matmul(OMP)  : ",real(elps_mtml_cpu),"(seconds)"
     write(*,*)'                                                               '
     write(*,*)'***************GPU Matrix Mul with cublasSgemm*****************'
     write(*,*)'                                                               '

     call start_timer_cpu("mtml_gpu")
     call gettimeofday_c(st_ZM_r)

     call my_SMatmul("N", "N", ipart, nrot, nk, &
                     alpha,                       &
                     pft1, lda,                   &
                     pft2, ldb,                   &
                     beta,                        &
                     prod_gpu,ldc                 )
     call gettimeofday_c(et_ZM_r)
     call elapsed_time_c(st_ZM_r,et_ZM_r,elps_mtml_gpu)
     call stop_timer_cpu("mtml_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_polarft_multi_gpus_gpu_c_ = -1"
     
     write(*,*)"Elapsed time cublasSgemm  : ",real(elps_mtml_gpu),"(seconds)"
     write(*,*)'                                                               '
     write(*,*)'*********Comparing product from matmul vs Sgemm (GPU)**********'
     write(*,*)'                                                               '

     write(*,'(10x,a,5x,a,5x,a,4x,a)') &
          "ipart","irot","(prod_cpu)","(prod_gpu)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
              write(*,*)i, irot, prod_cpu(i,irot), prod_gpu(i,irot)
           end do
        end do

     sum_diff_prod = sum(prod_cpu-prod_gpu)
     write(*,*) "Diff      sum(prod_cpu-prod_gpu): ",sum_diff_prod
     write(*,*) "N    prod_{cpu,gpu}(npart,npart): ",size_prod/1.e6,"(Millions)"
     write(*,*)'***************************************************************'
     speedup = elps_mtml_cpu / elps_mtml_gpu
     gflops = gflops / elps_mtml_gpu
     write(*,'(3x,a,2x,a,6x,a)') "Matmul CPU (secs)","Sgemm GPU (secs)","Gflp/s"
     write(*,'(x,f15.5,3x,f15.5,7x,f10.5)')elps_mtml_cpu,elps_mtml_gpu,gflops
     write(*,*)
     write(*,'(x,a,x,f10.3)')" Speed up matmul vs cublasSgemm: ",speedup
     write(*,*)'***************************************************************'

#if defined (BENCH)
     if( ibench_write ) then
        open(4,file='SMatmult_2D.asc',status='unknown',position='append')
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

  end do !end of do ipart = start_npart, npart, istep
  
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
  

  !shutting down the timers
  call stop_Alltimers_cpu()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_deviceQuery_gpu()

contains

end program testing_SMatmult_gpu
!//////////////////// Helper subroutine and functions //////////////////////////
!//
!*******************************************************************************
!    Subroutine to run sanity checks on the array of data structure passed GPU
!
!*******************************************************************************
!
subroutine Sanity_check_a_devD_gpu(a_devD,ndev)
  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_deviceQuery_gpu
  implicit none

  type(deviceDetails)             :: a_devD(*)
  integer :: ndev
  !local variables
  character,pointer,dimension(ndev) :: devname(:)
  integer :: idev

  !start of the execution commands
  
  write(*,*)"ndev                : ",a_devD(1:ndev)%ndev
  write(*,*)"d_ver               : ",a_devD(1:ndev)%d_ver
  write(*,*)"d_runver            : ",a_devD(1:ndev)%d_runver
  write(*,*)"tot_global_mem_MB   : ",a_devD(1:ndev)%tot_global_mem_MB
  write(*,*)"tot_global_mem_bytes: ",a_devD(1:ndev)%tot_global_mem_bytes
  write(*,*)"nMultiProc          : ",a_devD(1:ndev)%nMultiProc
  write(*,*)"ncc_per_mp          : ",a_devD(1:ndev)%ncc_per_mp
  write(*,*)"ncc                 : ",a_devD(1:ndev)%ncc
  write(*,*)"is_SMsuitable       : ",a_devD(1:ndev)%is_SMsuitable
  write(*,*)"nregisters_per_blk  : ",a_devD(1:ndev)%nregisters_per_blk
  write(*,*)"warpSze             : ",a_devD(1:ndev)%warpSze
  write(*,*)"maxthreads_per_mp   : ",a_devD(1:ndev)%maxthreads_per_mp
  write(*,*)"maxthreads_per_blk  : ",a_devD(1:ndev)%maxthreads_per_blk
  write(*,*)"is_ecc              : ",a_devD(1:ndev)%is_ecc
  write(*,*)"is_p2p              : ",a_devD(1:ndev)%is_p2p(1,0)

  return
end subroutine Sanity_check_a_devD_gpu
!*******************************************************************************
!    Subroutine to run sanity checks on the data structure passed GPU
!
!*******************************************************************************
!
subroutine Sanity_check_gpu(devQ, devD)
  use simple_defs
  use simple_cuda_defs
  use greeting_version
  use simple_deviceQuery_gpu
  implicit none

  type(deviceQuery_gpu)           :: devQ
  type(deviceDetails)             :: devD
  !local variables
  integer                         :: ndev

  !start of the execution commands

  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(devQ) and the data structure(devD) '
  write(*,*)'******************************************************************'
  write(*,*)'CUDA driver Ver.              (devQ):',devQ%get_cuda_drv_version()
  write(*,*)'CUDA driver Ver.              (devD):',devD%d_ver
  if(devD%d_ver/=devQ%get_cuda_drv_version())call devQ%get_warning_dataStruct_gpu()
  ndev = devQ%get_devCnt()
  write(*,*)'CUDA driver run Ver.          (devD):',devD%d_runver
  write(*,*)'DeviceCount                   (devQ):',ndev
  write(*,*)'DeviceCount                   (devD):',devD%ndev
  if ( ndev /= devD%ndev) call devQ%get_warning_dataStruct_gpu()
  write(*,*)'Device name                   (devD): ',devQ%get_devname(0)
  write(*,*)'Total global Mem in MB        (devD):',devD%tot_global_mem_MB
  write(*,*)'Total global Mem in bytes     (devD):',devD%tot_global_mem_bytes
  write(*,*)'nMultiProc                    (devD):',devD%nMultiProc
  write(*,*)'ncc_per_mp                    (devD):',devD%ncc_per_mp
  write(*,*)'ncc                           (devD):',devD%ncc
  write(*,*)'is_SMsuitable                 (devD):',devD%is_SMsuitable
  write(*,*)'nregisters_per_blk            (devD):',devD%nregisters_per_blk
  write(*,*)'warpSze                       (devD):',devD%warpSze
  write(*,*)'maxthreads_per_mp             (devD):',devD%maxthreads_per_mp
  write(*,*)'maxthreads_per_blk            (devD):',devD%maxthreads_per_blk
  write(*,*)'is_ecc                        (devD):',devD%is_ecc
  write(*,*)'******************************************************************'
  
  return
end subroutine Sanity_check_gpu
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
