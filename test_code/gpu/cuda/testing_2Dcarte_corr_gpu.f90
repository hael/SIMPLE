! ============================================================================
! Name        : testing_2Dcarte_corr_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 19th of April 2016
! Description : tests the cartesian corr calculation
!             :
! ============================================================================
!
program testing_2Dcarte_corr_gpu
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use simple_timing
  use simple_image,       only: image
  use simple_ft_expanded
  use simple_syscalls
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_deviceQuery_gpu
  use simple_systemQuery_cpu
  use simple_math, only: calc_corr, csq
  use simple_yaml_output
  use simple_yaml_strings

  !$ use omp_lib
  !$ use omp_lib_kinds

  implicit none
  logical, parameter :: ibench_local=.true.       !< benchmark result or not
  logical, parameter :: ibench_write_local=.true. !write benchmark result or not

  logical, parameter :: debug=.true.          !< debug indicator
  logical, parameter :: debug_cpu=.false.     !false/true
  logical, parameter :: debug_high=.true.     !true/false
  logical, parameter :: debug_write=.false.   !false/true
  logical, parameter :: debug_write_C=.false. !false/true

#define devptr_t integer*8

  type(polar_corr_calc)         :: s_carte
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD  
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
  type(image)                   :: img1, img2
  type(ft_expanded)             :: ftexp1, ftexp2

  integer, parameter            :: npart=2096
  integer, parameter            :: start_npart=2096, istep=1000

  integer                       :: ikernel = 6           !Kernel option
  integer                       :: threadsPerBlock = 512 !threadsPerBlock
  integer                       :: nx_3D=16, ny_3D=16, nz_3D=4 !3D mesh

  integer                       :: omp_nthr=1           !openMP variables
  integer, parameter            :: NITS=1

  !local variables
  type(deviceDetails),allocatable :: a_devD(:)
  !gpu gear
  integer                       :: vx, vy, vz !3D volume
  integer                       :: ndev
  integer                       :: size_cmat1, size_cmat2  
  integer                       :: rc
  integer                       :: lda
  integer                       :: ldb,ldc
  real(sp)                      :: alpha
  real(sp)                      :: r_gpu
  !CPU variables
  integer, allocatable          :: phys(:,:,:,:)
  logical, allocatable          :: lmsk(:,:,:)
  integer                       :: lims(3,2), inits
  real, parameter               :: lp=8.0
  real                          :: corr
  real(sp)                      :: r_cpu,sumasq_cpu,sumbsq_cpu
  !GPU variables
  complex(sp), allocatable      :: cmat1(:,:,:)       !< Fourier components
  complex(sp), allocatable      :: cmat2(:,:,:)       !< Fourier components

  !CUDA err variable for the return function calls
  integer                       :: err
  !timer variables
  double precision              :: elps_C
  double precision,dimension(2) :: st_C_r, et_C_r
  double precision              :: elps_L
  double precision,dimension(2) :: st_L_r, et_L_r

  double precision              :: elps_F
  double precision,dimension(2) :: st_F_r, et_F_r

  double precision              :: elps_O
  double precision,dimension(2) :: st_O_r, et_O_r
  double precision              :: elps_R
  double precision,dimension(2) :: st_R_r, et_R_r
  real(sp)                      :: speedup_OR
  real(sp)                      :: speedup_RF
  real(sp)                      :: speedup_OF
  real(sp)                      :: speedup_RC

  !indexers
  integer                       :: idev
  integer                       :: ivx,ivy,ivz
  !index variable major
  integer                       :: ipart

  !infinity handler
  integer :: inf
  real :: infinity
  equivalence (inf,infinity) !Stores two variable at the same address
  data inf/z'7f800000'/      !Hex for +Infinity

  !functions calls
  integer :: get_dev_count_c
  integer :: get_polarft_corr_gpu_c
  integer :: get_carte2d_ftext_corr_gpu_c

  !start of the execution commands
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

  !start of the execution commands
  call timestamp()
  call start_Alltimers_cpu()
  !$ call omp_set_num_threads(omp_nthr)

  !starting the cuda environment
  call simple_cuda_init(err,devQ,devD)
  if (err .ne. RC_SUCCESS ) write(*,*) 'cublas init failed'

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   System fills in the object(sysQ) and the data structure(hstD)  '
  write(*,*)'******************************************************************'
  call sysQ%new_systemQuery_cpu(hstD)
  !call Sanity_check_cpu(sysQ, hstD)
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

  !call Sanity_check_a_devD_gpu(a_devD,ndev)

  do ipart = start_npart, npart, istep

!*******************************************************************************
!     now testign the corr function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************
     write(*,*)'                                                               '
     write(*,*)'***************CPU corr****************************************'
     write(*,*)'                                                               '

     ! make random images
     call img1%new([ipart,ipart,1], 1.77)
     call img1%ran
     call img2%new([ipart,ipart,1], 1.77)
     call img2%ran
     ! prepare objects for corrcalc
     call ftexp1%new(img1, lp)
     call ftexp2%new(img2, lp)
     lims = img1%loop_lims(1,lp)
     ! time the different routines
     call start_timer_cpu("old")
     call gettimeofday_c(st_O_r)
     do inits=1,NITS
        corr  = img1%corr(img2, lp)
     end do
     call gettimeofday_c(et_O_r)
     call elapsed_time_c(st_O_r,et_O_r,elps_O)
     call stop_timer_cpu("old")

     call start_timer_cpu("recast")
     call gettimeofday_c(st_R_r)
     do inits=1,NITS
        corr  = ftexp1%corr(ftexp2)
     end do
     call gettimeofday_c(et_R_r)
     call elapsed_time_c(st_R_r,et_R_r,elps_R)
     call stop_timer_cpu("recast")

     write(*,*)'                                                               '
     write(*,*)'************Allocating and initialising data for GPU***********'
     write(*,*)'                                                               '
     
     lda = ipart
     ldb = ipart
     ldc = ipart
!*******************************************************************************
!    Data allocation and initialisation of the matrices for testing purposes.
!
!*******************************************************************************
!     allocate(cmat1(lims(1,1):lims(1,2), lims(2,1):lims(2,2), &
!                    lims(3,1):lims(3,2))                      )
!     allocate(cmat2(lims(1,1):lims(1,2), lims(2,1):lims(2,2), &
!                    lims(3,1):lims(3,2))                      )

     write(*,*)"size of cmat1: ",size(cmat1,kind=sp)
     write(*,*)"size of cmat2: ",size(cmat2,kind=sp)
     write(*,*)"lims: ", lims(1,1), lims(1,2), &
          lims(2,1) , lims(2,2), &
          lims(3,1) , lims(3,2)

     vx = lims(1,2) - lims(1,1) + 1
     vy = lims(2,2) - lims(2,1) + 1
     vz = lims(3,2) - lims(3,1) + 1

     write(*,*) vx,vy,vz
     
     allocate(cmat1(vx, vy, vz))
     allocate(cmat2(vx, vy, vz))

     !initialisatino of the cmat data
     cmat1 = 0.0
     cmat2 = 0.0
     
!*******************************************************************************
!    Filling and printing cmat array with random numbers using test functions
!
!*******************************************************************************
     call getCRandomGaussianDistr_3D(vx, vy, vz, cmat1)
     call getCRandomGaussianDistr_3D(vx, vy, vz, cmat2)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a)')"ivx","ivy","ivz", &
                                              "(cmat1)","(cmat2)"
     do ivx=1,vx - (vx-2)
        do ivy=1,vy - (vy-2)
           do ivz=1,vz - (vz-1)
              write(*,*)ivx,ivy,ivz,cmat1(ivx,ivy,ivz),cmat2(ivx,ivy,ivz)
           end do
        end do
     end do
     write(*,'(50x,a,12x,a)')"conjg(cmat2)","real(cmat1(ivx,ivy,ivz)*conjg(cmat2(ivx,ivy,ivz)))"
     do ivx=1,vx - (vx-2)
        do ivy=1,vy - (vy-2)
           do ivz=1,vz - (vz-1)
              write(*,*)ivx, ivy, ivz, conjg(cmat2(ivx,ivy,ivz)) , &
                   real( cmat1(ivx,ivy,ivz)*conjg(cmat2(ivx,ivy,ivz)) ), &
                   real(cmat1(ivx,ivy,ivz)) * real(cmat2(ivx,ivy,ivz)) + &
                   imag(cmat1(ivx,ivy,ivz)) * imag(cmat2(ivx,ivy,ivz)) 
           end do
        end do
     end do
!*******************************************************************************
!     now testign the corr function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************
     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     Now testing the corr function                             '
     write(*,*)'***************************************************************'
     size_cmat1 = vx*vy*vz
     size_cmat2 = vx*vy*vz
     write(*,*)
     write(*,'(3x,a,3x,a,6x,a,6x,a,6x,a,6x,a)') &
          "npart","ipart","vx","vy","vz","OpenMP thrd"
     write(*,'(6(3x,i5))') npart,  ipart, vx, vy, vz, omp_nthr
     write(*,*) "N cmat1(npart,nrot/2,nk): ",size_cmat1
     write(*,*) "N cmat1(npart,nrot/2,nk): ",size_cmat1/1.e6,"(Millions)"
     write(*,*) "N cmat2(nx,ny,nz): ",size_cmat2
     write(*,*) "N cmat2(nx,ny,nz): ",size_cmat2/1.e6,"(Millions)"
     write(*,*) " Ratio of size (CMAT2 and CMAT1): ",size_cmat2/real(size_cmat1)

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     alpha = 1.0
     s_carte%r_polar         = 1.342
     s_carte%sumasq_polar    = 2.132
     s_carte%sumbsq_polar    = 3.123
     s_carte%ikrnl           = ikernel
     s_carte%threadsPerBlock = threadsPerBlock
     s_carte%nx              = nx_3D
     s_carte%ny              = ny_3D
     s_carte%nz              = nz_3D

     write(*,*)'                           F    N                              '
     r_gpu = 1.0
     call start_timer_cpu("recastFN_gpu")
     call gettimeofday_c(st_F_r)
     err = get_polarft_corr_gpu_c(a_devD,                &
                                  s_carte,"F","N",       &
                                  r_gpu,                 &
                                  cmat1,cmat2,           &
                                  vx,vy,vz,              &
                                  lda,ldb,ldc,alpha,     &
                                  s_bench, s_debug_gpu)
     call gettimeofday_c(et_F_r)
     call elapsed_time_c(st_F_r,et_F_r,elps_F)
     call stop_timer_cpu("recastFN_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_polarft_corr_gpu_c=",err

     r_gpu = calc_corr(s_carte%r_polar,s_carte%sumasq_polar*s_carte%sumbsq_polar)
     write(*,*)'                                                               '
     write(*,*)'*************Comparing correlator from CPU vs GPU**************'
     write(*,*)'                                                               '

     call start_timer_cpu("vs")
     call gettimeofday_c(st_L_r)
     ! corr is real part of the complex mult btw 1 and 2*
     r_cpu = 0.0
     r_cpu = sum(real(cmat1*conjg(cmat2)))
     write(*,*)"Hadamard r_cpu:",r_cpu,     "     r_gpu:",s_carte%r_polar
     ! normalisation terms
     sumasq_cpu = sum(csq(cmat1))
     sumbsq_cpu = sum(csq(cmat2))
     write(*,*)"sum sumasq_cpu:",sumasq_cpu,"sumasq_gpu:",s_carte%sumasq_polar
     write(*,*)"sum sumbsq_cpu:",sumbsq_cpu,"sumbsq_gpu:",s_carte%sumbsq_polar
     ! finalise the correlation coefficient
     r_cpu = calc_corr(r_cpu,sumasq_cpu*sumbsq_cpu)
     write(*,*)"after calc_corr r_cpu : ",r_cpu," r_gpu : ",r_gpu
     call gettimeofday_c(et_L_r)
     call elapsed_time_c(st_L_r,et_L_r,elps_L)
     call stop_timer_cpu("vs")
    
     speedup_OR = elps_O/elps_R
     speedup_RF = elps_R/elps_F
     speedup_OF = elps_O/elps_F
     write(*,*)'***************************************************************'
     write(*,'(x,a,12x,a,11x,a,7x,a)')"Correlations","old","re-cast","re-cast(GPU)"
     write(*,*)"---------------------------------------------------------------"
     write(*,'(16x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8)') elps_O,elps_R,elps_F,elps_L
     write(*,*)'***************************************************************'
     write(*,'(3x,a,3x,a,6x,a,6x,a,6x,a,6x,a)') &
          "npart","ipart","vx","vy","vz","OpenMP thrd"
     write(*,'(6(3x,i5))') npart,  ipart, vx, vy, vz, omp_nthr
     write(*,*) " Speed up from old to recast    : ",speedup_OR
     if (speedup_OR<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speed up from recast to GPU(F) : ",speedup_RF
     if (speedup_RF<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speed up from old to GPU(F)    : ",speedup_OF
     if (speedup_OF<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*)'***************************************************************'

#if defined (BENCH)
     if( ibench_write ) then
        open(4,file='2Dcarte_corr_Al.asc',status='unknown',position='append')
        write(4,'(4(1x,i5),6(2x,f15.8))')ipart,vx,vy,vz,elps_O,elps_R,elps_F, &
             speedup_OR, speedup_RF, speedup_OF
        close(4)
     end if
#endif
!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************

   deallocate(cmat1)
   deallocate(cmat2)

  end do

  !freeing ressources on host for the data structure 
  deallocate(a_devD)

!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
  !shutting down the environment
  call simple_cuda_shutdown()
  !shutting down the timers
  call stop_Alltimers_cpu()
  
end program testing_2Dcarte_corr_gpu
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
