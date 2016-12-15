! ============================================================================
! Name        : testing_2DshCar_corr_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 19th of April 2016
! Description : tests the cartesian corr calculation
!             :
! ============================================================================
!
program testing_2DshCar_corr_gpu
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
  logical, parameter :: ibench_local=.true.      !< benchmark result or not
  logical, parameter :: ibench_write_local=.true.!write benchmark result or not

  logical, parameter :: debug=.false.          !< debug indicator
  logical, parameter :: debug_cpu=.false.     !false/true
  logical, parameter :: debug_high=.false.     !true/false
  logical, parameter :: debug_write=.false.   !false/true
  logical, parameter :: debug_write_C=.false. !false/true

#define devptr_t integer*8

  type(polar_corr_calc)         :: s_carte,s_carteCN,s_carteCF
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD  
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
  type(image)                   :: img1, img2
  type(ft_expanded)             :: ftexp1, ftexp2

  integer, parameter            :: npart=30096
  integer, parameter            :: start_npart=1096, istep=1000

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
  real, parameter               :: lp=8.0, shvec(3)=[0.,0.,0.]
  real                          :: corr
  !cpu
  real(sp)                      :: r_cpu, sumasq_cpu, sumbsq_cpu
  real(sp)                      :: rs_cpu, s_sumasq_cpu, s_sumbsq_cpu
  !gpu
  real(sp)                      :: rs_gpu, s_sumasq_gpu, s_sumbsq_gpu
  !CPU and GPU variables
  complex(sp), allocatable      ::   cmat1(:,:,:) !< Fourier components
  complex(sp), allocatable      ::   cmat2(:,:,:) !< Fourier components
  complex(sp), allocatable      ::   shmat(:,:,:)
  complex(sp), allocatable      :: cmat2sh(:,:,:)

  !CUDA err variable for the return function calls
  integer                       :: err

  !timer variables
  double precision              :: elps_CF         !rescast-shifted GPU (CF)
  double precision,dimension(2) :: st_CF_r, et_CF_r

  double precision              :: elps_CN         !rescast-shifted GPU (CN)
  double precision,dimension(2) :: st_CN_r, et_CN_r

  double precision              :: elps_L         !recast local
  double precision,dimension(2) :: st_L_r, et_L_r

  double precision              :: elps_F         !rescast GPU (FN)
  double precision,dimension(2) :: st_F_r, et_F_r

  double precision              :: elps_O         !old ft_expanded
  double precision,dimension(2) :: st_O_r, et_O_r

  double precision              :: elps_R         !recast ft_expanded
  double precision,dimension(2) :: st_R_r, et_R_r

  double precision              :: elps_S         !rescast-shifted ft_expanded
  double precision,dimension(2) :: st_S_r, et_S_r

  double precision              :: elps_H         !recast-shifted local
  double precision,dimension(2) :: st_H_r, et_H_r
  !Relative speed up for the diffenrent methods 
  real(sp)                      :: speedup_OR
  real(sp)                      :: speedup_RF
  real(sp)                      :: speedup_OF
  real(sp)                      :: speedup_RCN
  real(sp)                      :: speedup_SCN
  real(sp)                      :: speedup_HCN
  real(sp)                      :: speedup_RCF
  real(sp)                      :: speedup_SCF
  real(sp)                      :: speedup_HCF
  real(sp)                      :: speedup_CFN

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

!*******************************************************************************
!     now testign the corr function on both CPU(OpenMP) vs GPU single thread
!
!*******************************************************************************

#if defined (BENCH)
   if( ibench_write ) then
      open(4,file='2Dcarte_corr_Al.asc',status='unknown',position='append')
      write(4,'(1x,a,7x,a,34x,a,45x,a)')"[x,x,1]","lims(i,j)",&
           "<-----------Correlations(secs)(old--->(lc))------>",&
           "<------------speedup relators(OR--->CFN)--------->"
      write(4,'(4x,a,5x,a,4x,a,4x,a,6x,a,6x,a,4x,a,4x,a,2x,a,1x,a,3x,a,3x,a,7x,a,9x,a,9x,a,7(8x,a))')&
           "x","vx","vy","vz",                  &
           "Old","Rcst_ft","ft-shft",           &
           "GPU(FN)","shft-GPU(CF)","shft(CN)", &
           "CPU(lc)","shft(lc)",                &
           "OR","RF","OF",                      &
           "RCF","SCF","HCF",                   &
           "RCN","SCN","HCN",                   &
           "CFN"
      close(4)
   end if
#endif
  
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

     call start_timer_cpu("shifted")
     call gettimeofday_c(st_S_r)
     do inits=1,NITS
        corr = ftexp1%corr_shifted(ftexp2, shvec)
     end do
     call gettimeofday_c(et_S_r)
     call elapsed_time_c(st_S_r,et_S_r,elps_S)
     call stop_timer_cpu("shifted")
     
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
     allocate(shmat(vx, vy, vz))
     allocate(cmat2sh(vx, vy, vz))

     !initialisatino of the cmat data
     cmat1 = 0.0
     cmat2 = 0.0
     shmat = 0.0
     cmat2sh = 0.0
     
!*******************************************************************************
!    Filling and printing cmat array with random numbers using test functions
!
!*******************************************************************************
     call getCRandomGaussianDistr_3D(vx, vy, vz, cmat1)
     call getCRandomGaussianDistr_3D(vx, vy, vz, cmat2)
     call getCRandomGaussianDistr_3D(vx, vy, vz, shmat)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a,30x,a)')"ivx","ivy","ivz", &
                                              "(cmat1)","(cmat2)","(shmat)"
     do ivx=1,vx - (vx-2)
        do ivy=1,vy - (vy-2)
           do ivz=1,vz - (vz-1)
              write(*,*)ivx,ivy,ivz,cmat1(ivx,ivy,ivz),cmat2(ivx,ivy,ivz),shmat(ivx,ivy,ivz)
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
     call start_timer_cpu("u_recastFN_gpu")
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
     call stop_timer_cpu("u_recastFN_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_polarft_corr_gpu_c=",err

     r_gpu = calc_corr(s_carte%r_polar,s_carte%sumasq_polar*s_carte%sumbsq_polar)

     write(*,*)'                           C    F                              '

     s_carteCF%r_polar         = 1.342
     s_carteCF%sumasq_polar    = 2.132
     s_carteCF%sumbsq_polar    = 3.123
     s_carteCF%ikrnl           = ikernel
     s_carteCF%threadsPerBlock = threadsPerBlock
     s_carteCF%nx              = nx_3D
     s_carteCF%ny              = ny_3D
     s_carteCF%nz              = nz_3D

     rs_gpu = 1.0

     call start_timer_cpu("v_recastCF_gpu")
     call gettimeofday_c(st_CF_r)
     err = get_carte2d_ftext_corr_gpu_c(a_devD,                     &
                                        s_carteCF,"C","F",          &
                                        rs_gpu, shmat,              &        
                                        cmat1,cmat2,                &
                                        vx,vy,vz,                   &
                                        lda,ldb,ldc,alpha,          &
                                        s_bench, s_debug_gpu) 
     call gettimeofday_c(et_CF_r)
     call elapsed_time_c(st_CF_r,et_CF_r,elps_CF)
     call stop_timer_cpu("v_recastCF_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_carte2d_ftext_corr_gpu_c=",err

     rs_gpu = calc_corr(s_carteCF%r_polar,s_carteCF%sumasq_polar*s_carteCF%sumbsq_polar)

     write(*,*)'                           C    N                              '

     s_carteCN%r_polar         = 1.342
     s_carteCN%sumasq_polar    = 2.132
     s_carteCN%sumbsq_polar    = 3.123
     s_carteCN%ikrnl           = ikernel
     s_carteCN%threadsPerBlock = threadsPerBlock
     s_carteCN%nx              = nx_3D
     s_carteCN%ny              = ny_3D
     s_carteCN%nz              = nz_3D

     rs_gpu = 1.0

     call start_timer_cpu("v_recastCN_gpu")
     call gettimeofday_c(st_CN_r)
     err = get_carte2d_ftext_corr_gpu_c(a_devD,                     &
                                        s_carteCN,"C","N",          &
                                        rs_gpu, shmat,              &        
                                        cmat1,cmat2,                &
                                        vx,vy,vz,                   &
                                        lda,ldb,ldc,alpha,          &
                                        s_bench, s_debug_gpu) 
     call gettimeofday_c(et_CN_r)
     call elapsed_time_c(st_CN_r,et_CN_r,elps_CN)
     call stop_timer_cpu("v_recastCN_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_carte2d_ftext_corr_gpu_c=",err

     rs_gpu = calc_corr(s_carteCN%r_polar,s_carteCN%sumasq_polar*s_carteCN%sumbsq_polar)

     write(*,*)'                                                               '
     write(*,*)'*************Comparing correlator from CPU vs GPU**************'
     write(*,*)'                                                               '

     write(*,*)'***************CPU corr vs GPU(FN)*****************************'
     call start_timer_cpu("vs_corr")
     call gettimeofday_c(st_L_r)
     ! corr is real part of the complex mult btw 1 and 2*
     r_cpu = 0.0
     r_cpu = sum(real(cmat1*conjg(cmat2)))
     write(*,*)"Hadr r_cpu:",r_cpu,     "     r_gpu:",s_carte%r_polar
     ! normalisation terms
     sumasq_cpu = sum(csq(cmat1))
     sumbsq_cpu = sum(csq(cmat2))
     write(*,*)"sumasq_cpu:",sumasq_cpu,"sumasq_gpu:",s_carte%sumasq_polar
     write(*,*)"sumbsq_cpu:",sumbsq_cpu,"sumbsq_gpu:",s_carte%sumbsq_polar
     ! finalise the correlation coefficient
     r_cpu = calc_corr(r_cpu,sumasq_cpu*sumbsq_cpu)
     write(*,*)"after calc_corr r_cpu : ",r_cpu," r_gpu : ",r_gpu
     call gettimeofday_c(et_L_r)
     call elapsed_time_c(st_L_r,et_L_r,elps_L)
     call stop_timer_cpu("vs_corr")

     write(*,*)'***************CPU corrshifted vs GPU(CN)**********************'

     rs_cpu = 0.0
     s_sumasq_cpu = 0.0
     s_sumbsq_cpu = 0.0

     call start_timer_cpu("vs_corr_shifted")
     call gettimeofday_c(st_H_r)
     cmat2sh = cmat2*shmat
     do ivx=1,vx - (vx-2)
        do ivy=1,vy - (vy-2)
           do ivz=1,vz - (vz-1)
              write(*,*)ivx,ivy,ivz,cmat2sh(ivx,ivy,ivz)
           end do
        end do
     end do
     ! corr is real part of the complex mult btw 1 and 2*
     rs_cpu = sum(real(cmat1*conjg(cmat2sh)))
     write(*,*)"Hadmr rs_cpu:",rs_cpu,     "    rs_gpu:",s_carteCN%r_polar
     ! normalisation terms
     s_sumasq_cpu = sum(csq(cmat1))
     s_sumbsq_cpu = sum(csq(cmat2sh))
     ! finalise the correlation coefficient
     rs_cpu = calc_corr(rs_cpu,s_sumasq_cpu*s_sumbsq_cpu)
     write(*,*)"s_sumasq_cpu:",s_sumasq_cpu,"sumasq_gpu:",s_carteCN%sumasq_polar
     write(*,*)"s_sumbsq_cpu:",s_sumbsq_cpu,"sumbsq_gpu:",s_carteCN%sumbsq_polar
     write(*,*)"after calc_corr rs_cpu : ",rs_cpu," rs_gpu : ",rs_gpu
     call gettimeofday_c(et_H_r)
     call elapsed_time_c(st_H_r,et_H_r,elps_H)
     call stop_timer_cpu("vs_corr_shifted")
     
     speedup_OR = elps_O/elps_R
     speedup_RF = elps_R/elps_F
     speedup_OF = elps_O/elps_F

     speedup_RCN = elps_R/elps_CN
     speedup_SCN = elps_S/elps_CN
     speedup_HCN = elps_H/elps_CN

     speedup_RCF = elps_R/elps_CF
     speedup_SCF = elps_S/elps_CF
     speedup_HCF = elps_H/elps_CF

     speedup_CFN = elps_CF/elps_CN

     write(*,*)'***************************************************************'
     write(*,'(x,a,5x,a,11x,a,6x,a,2x,a,2x,a,2x,a)')&
          "Corrs:","old","re-cast","rcst(GPU)(FN)",&
          "rcstsh(GPU)(CF)","rcstsh(GPU)(CN)","rcstsh(CPU)(ft)",&
          "recastsh(CPU)"
     write(*,*)"---------------------------------------------------------------"
     write(*,'(3x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8,1x,f15.8)') &
          elps_O,elps_R,elps_F, &
          elps_CF,elps_CN,elps_S,elps_H
     write(*,*)'***************************************************************'
     write(*,'(3x,a,3x,a,6x,a,6x,a,6x,a,6x,a)') &
          "npart","ipart","vx","vy","vz","OpenMP thrd"
     write(*,'(6(3x,i5))') npart, ipart, vx, vy, vz, omp_nthr

     write(*,*) " Speed up from old to recast      : ",speedup_OR
     if (speedup_OR<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speed up from recast to GPU(F)   : ",speedup_RF
     if (speedup_RF<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speed up from old to GPU(FN)     : ",speedup_OF
     if (speedup_OF<1.0) write(*,*)"speedup < 1, try bigger Volume"

     write(*,*) " Speed up recast to GPU(CF)       : ",speedup_RCF
     if (speedup_RCF<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speedup recast-shft_ft to GPU(CF): ",speedup_SCF
     if (speedup_SCF<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speedup recast-shft to GPU(CF)   : ",speedup_HCF
     if (speedup_HCF<1.0) write(*,*)"speedup < 1, try bigger Volume"

     write(*,*) " Speed up recast to GPU(CN)       : ",speedup_RCN
     if (speedup_RCN<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speedup recast-shft_ft to GPU(CN): ",speedup_SCN
     if (speedup_SCN<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*) " Speedup recast-shft to GPU(CN)   : ",speedup_HCN
     if (speedup_HCN<1.0) write(*,*)"speedup < 1, try bigger Volume"

     write(*,*) " Speedup GPU(CF) to GPU(CN)       : ",speedup_CFN
     if (speedup_CFN<1.0) write(*,*)"speedup < 1, try bigger Volume"
     write(*,*)'***************************************************************'

#if defined (BENCH)
     if( ibench_write ) then
        open(4,file='2Dcarte_corr_Al.asc',status='unknown',position='append')
        write(4,'(4(1x,i5),18(1x,f10.5))')         &
             ipart,vx,vy,vz,                       &
             elps_O, elps_R, elps_S, elps_F,elps_CF,elps_CN, elps_L, elps_H, &
             speedup_OR, speedup_RF, speedup_OF,   &
             speedup_RCF,speedup_SCF, speedup_HCF, &
             speedup_RCN,speedup_SCN, speedup_HCN, &
             speedup_CFN
        close(4)
     end if
#endif
!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************

   deallocate(cmat1)
   deallocate(cmat2)
   deallocate(shmat)
   deallocate(cmat2sh)

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
  
end program testing_2DshCar_corr_gpu
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
