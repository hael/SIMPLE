! ============================================================================
! Name        : testing_MultiGPU_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 07th of March 2016
! Description : tests the Multi GPU environment within CUDA
!             :
! ============================================================================
!
program testing_MultiGPU_gpu
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
  use simple_random

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

  type(polar_corr_calc)         :: s_polar
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD  
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
!  type(cuda)                    :: cudaQ
  integer, parameter            :: npart=16
  integer, parameter            :: start_npart=16, istep=64
  integer, parameter            :: nrot = 400
  integer, parameter            :: nk = 100

  integer                       :: ikernel = 0           !Kernel option
  integer                       :: threadsPerBlock = 256 !threadsPerBlock
  integer                       :: nx_3D=16, ny_3D=16, nz_3D=4 !3D mesh
  
  integer                       :: omp_nthr=8           !openMP variables

  !local variables
  type(deviceDetails),allocatable :: a_devD(:)
  character(len=1),allocatable   :: devname(:,:)
  !gpu gear
  integer                       :: ndev
  integer                       :: rc

  integer                       :: refsz,ptclsz
  integer                       :: size_pft1, size_pft2  
  integer                       :: size_crmt3d, size_crmt2d
  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc
  real(sp)                      :: alpha

  !comparators
  real(sp)                      :: sum_diff_crmt2d
  real(sp)                      :: sum_diff_crmt3d
  !variables for the function 

  real(sp), allocatable         ::      sumZ_vec(:)
  real(sp), allocatable         ::   sqsums_refs(:)    !self%sqsums_refs(self%nrefs)
  real(sp), allocatable         ::  sqsums_ptcls(:)    !self%sqsums_ptcls(self%nptcls)
  real(sp), allocatable         ::     corrmat2d(:,:)  !cor matrix 2d final output
  integer, allocatable          ::     inplmat2d(:,:)  !inplane matrix 2d final output
  real(sp), allocatable         ::     corrmat3d(:,:,:)!cor matrix 3d final output
  real(sp), allocatable         :: hadamard_prod(:,:,:)!(self%nptcls,self%refsz,self%nk)
  complex(sp), allocatable      ::          pft1(:,:,:)
  complex(sp), allocatable      ::          pft2(:,:,:)
  
  !gpu variables for the calculation of the final output

  real(sp), allocatable         :: corrmat3d_gpu(:,:,:) !cor matrix 3d
  real(sp), allocatable         :: corrmat2d_gpu(:,:)   !cor matrix 2d
  integer, allocatable          :: inplmat2d_gpu(:,:)   !inplane matrix 2d 

  ! dimension arrays
  integer                       :: d1lim(2), d2lim(2)

  !timer variables

  real(sp)                      :: speedup

  double precision              :: elps_corr_ZM_gpu
  double precision,dimension(2) :: st_ZM_r, et_ZM_r

  double precision              :: elps_corr_cpu
  double precision,dimension(2) :: st_corr_cpu, et_corr_cpu

  double precision              :: elps_mze_ref_cpu
  double precision,dimension(2) :: st_mze_ref_cpu, et_mze_ref_cpu
  
  double precision              :: elps_mze_ptc_cpu
  double precision,dimension(2) :: st_mze_ptc_cpu, et_mze_ptc_cpu

  !indexers
  integer                       :: i,j,jpart,ik,irot,iptcl
  integer                       :: idev
  !index variable major
  integer                       :: ipart

  !infinity handler
  integer :: inf
  real :: infinity
  equivalence (inf,infinity) !Stores two variable at the same address
  data inf/z'7f800000'/      !Hex for +Infinity

  !functions calls
  integer :: get_dev_count_c
  integer :: get_polarft_multi_gpus_gpu_c
  integer :: findCudaDevice_c
  integer :: setDevice_c
  integer :: strlen

  !start of the execution commands

  !initialising the fixed seed for the random number generator
  call init_fixed_random_seed(1234)

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

  !$ call omp_set_num_threads(omp_nthr)

  !starting the cuda environment
  call simple_cuda_init(err,devQ,devD)
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
  allocate(devname(0:devD%ndev-1,0:DEVNAME_STRING_LENGTH))
  do idev = 0, ndev-1
     call devQ%new_deviceQuery_gpu(devD,idev)
     devname(idev,:) = devQ%get_devname(0)
     devname(idev,:) = devname(idev,0:strlen(devname(idev,:)))
     !call Sanity_check_gpu(devQ, devD)
     !mapping the data structures into an array
     a_devD(idev) = devD
  end do

  call Sanity_check_a_devD_gpu(a_devD,ndev)

  do ipart = start_npart, npart, istep

     lda = ipart
     ldb = nrot
     ldc = ipart

!*******************************************************************************
!    Data allocation and initialisation of the matrices for testing purposes.
!
!*******************************************************************************

     refsz  = nrot / 2    !size of reference (nrots/2)
     ptclsz = nrot*2      !size of particle (2*nrots)

     !cpu allocation
     allocate(     sumZ_vec(ipart              ))
     allocate(  sqsums_refs(ipart              ))
     allocate( sqsums_ptcls(ipart*2            ))
     allocate(         pft1(ipart,    refsz, nk)) !pfts_refs
     allocate(         pft2(ipart*2, ptclsz, nk)) !pfts_ptcls

     allocate(    corrmat2d(ipart, ipart       ))
     allocate(    inplmat2d(ipart, ipart       ))
     allocate(    corrmat3d(ipart, ipart, nrot ))
     allocate(hadamard_prod(ipart, refsz, nk   ))
     
     sumZ_vec          = 0.0
     pft1              = 0.0
     pft2              = 0.0
     sqsums_refs       = 0.0
     sqsums_ptcls      = 0.0
     corrmat3d         = 0.0
     hadamard_prod     = 0.0
     corrmat2d         = 0.0
     inplmat2d         = 0
     !gpu allocation
     allocate(corrmat2d_gpu(ipart,ipart       ))
     allocate(inplmat2d_gpu(ipart,ipart       ))
     allocate(corrmat3d_gpu(ipart, ipart, nrot))

     corrmat3d_gpu     = 0.0
     corrmat2d_gpu     = 0.0
     inplmat2d_gpu     = 0

!*******************************************************************************
!    Filling and printing PFTs array with random numbers using test functions
!
!*******************************************************************************
     call getCRandomGaussianDistr_3D( ipart,  refsz, nk, pft1)
     call getCRandomGaussianDistr_3D(2*ipart, ptclsz, nk, pft2)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a)')"ipart","irot","ik", &
                                              "(pft1)","(pft2)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot,ik,pft1(i,irot,ik), pft2(i,irot,ik)
           end do
        end do
     end do
     write(*,'(50x,a,12x,a)')"conjg(pft2)","real(pft1(i1,irot,ik)*conjg(pft2(i2,irot,ik)))"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot, ik, conjg(pft2(i,irot,ik)) , &
                   real( pft1(i,irot,ik) * conjg(pft2(i,irot,ik)) ), &
                   real(pft1(i,irot,ik)) * real(pft2(i,irot,ik))   + &
                   imag(pft1(i,irot,ik)) * imag(pft2(i,irot,ik)) 
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
     size_pft1 = ipart*refsz*nk
     size_pft2 = (2*ipart)*ptclsz*nk
     size_crmt2d = ipart * npart
     size_crmt3d = ipart * npart * nrot
     write(*,*)
     write(*,'(3x,a,3x,a,4x,a,6x,a,3x,a,3x,a,2x,a,2x,a)') &
          "npart","ipart","nrot","nk","refsz","ptclsz","(2*ipart)","OpenMP thrd"
     write(*,'(8(3x,i5))') npart,  ipart,  nrot,  nk,  refsz,  ptclsz,  (2*ipart), omp_nthr
     write(*,*) "N (Ref)   pft1(npart,nrot/2,nk): ",size_pft1
     write(*,*) "N (Ref)   pft1(npart,nrot/2,nk): ",size_pft1/1.e6,"(Millions)"
     write(*,*) "N (ptc) pft2(npart*2,nrot*2,nk): ",size_pft2
     write(*,*) "N (ptc) pft2(npart*2,nrot*2,nk): ",size_pft2/1.e6,"(Millions)"
     write(*,*) " Ratio of size  (PFT2 and PFT1): ",size_pft2 / real(size_pft1)
     write(*,*)'                                                               '
     write(*,*)'***************CPU corr****************************************'
     write(*,*)'                                                               '

     call gettimeofday_c(st_mze_ref_cpu)
     call memoize_sqsum_ref(sqsums_refs, pft1, ipart, refsz, nk)
     call gettimeofday_c(et_mze_ref_cpu)
     call elapsed_time_c(st_mze_ref_cpu,et_mze_ref_cpu,elps_mze_ref_cpu)

     call gettimeofday_c(st_mze_ptc_cpu)
     call memoize_sqsum_ptcl(sqsums_ptcls, pft2, (2*ipart), refsz, nk)
     call gettimeofday_c(et_mze_ptc_cpu)
     call elapsed_time_c(st_mze_ptc_cpu,et_mze_ptc_cpu,elps_mze_ptc_cpu)
     
     call start_timer_cpu("corr_cpu")
     call gettimeofday_c(st_corr_cpu)
     call get_polarft_gencorrall_cpu(corrmat3d, &
                                     hadamard_prod,                       &
                                     d1lim, d2lim, s_polar, sumZ_vec,     &
                                     pft1,pft2,                           &
                                     sqsums_refs,sqsums_ptcls,            &
                                     ipart,nrot,nk,                       &
                                     lda,ldb,ldc,alpha)


     !$omp parallel default(shared) private(iptcl)
     !$omp workshare
     corrmat2d = maxval(corrmat3d, dim=3)
     inplmat2d = maxloc(corrmat3d, dim=3)
     !$omp end workshare
     !$omp end parallel

     call gettimeofday_c(et_corr_cpu)
     call elapsed_time_c(st_corr_cpu,et_corr_cpu,elps_corr_cpu)
     call stop_timer_cpu("corr_cpu")

     where ( corrmat2d == infinity ) corrmat2d = 0.0
     where ( corrmat3d ==  infinity ) corrmat3d = 0.0
     where ( corrmat3d == -infinity ) corrmat3d = 0.0

     write(*,*)"E[sqsums_refs(:)] is    : ",sum(sqsums_refs(:))/ipart
     write(*,'(x,a,i4,x,a,i4,x,a)') &
               "Over indx from d1lim(:) : [",d1lim(1),": ",d1lim(2),"]"
     write(*,*)"E[sqsums_ptcls(:)] is   : ",sum(sqsums_ptcls(d1lim(1):d1lim(2)))/ipart
     write(*,*)"Elapsed time    corr_cpu: ",real(elps_corr_cpu),   "(seconds)"
     write(*,*)"Elapsed time mze_ref_cpu: ",real(elps_mze_ref_cpu),"(seconds)"
     write(*,*)"Elapsed time mze_ptc_cpu: ",real(elps_mze_ptc_cpu),"(seconds)"

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     alpha = 1.0
     s_polar%r_polar         = 1.342
     s_polar%sumasq_polar    = 2.132
     s_polar%sumbsq_polar    = 3.123
     s_polar%ikrnl           = ikernel
     s_polar%threadsPerBlock = threadsPerBlock
     s_polar%nx              = nx_3D
     s_polar%ny              = ny_3D
     s_polar%nz              = nz_3D
     
     write(*,*)'                           Z    M                              '
     write(*,*) "looking for the best device to compute on..."
     idev = findCudaDevice_c(idev)
     write(*,*)"Device[",idev,"]:",devname(idev,:)," with ",a_devD(idev)%ncc," cores"
     err = setDevice_c(idev)
     write(*,'(3x,a,3x,a,4x,a,6x,a,3x,a,3x,a,2x,a,2x,a,2x,a,2x,a,2x,a)') &
          "npart","ipart","nrot","nk","refsz","ptclsz","(2*ipart)",&
          "nx, ny, nz","ThreadsPerBlock","Kernel"
     write(*,'(7(3x,i5),5x,i3,x,i3,x,i3,8x,i3,8x,i2)') &
          npart, ipart, nrot, nk, refsz, ptclsz, (2*ipart), &
          s_polar%nx,s_polar%ny,s_polar%nz, &
          s_polar%threadsPerBlock, &
          s_polar%ikrnl
     call start_timer_cpu("corr_ZM_gpu")
     call gettimeofday_c(st_ZM_r)
     err = get_polarft_multi_gpus_gpu_c(a_devD,                     &
                                        s_polar,"Z","M",            &
                                        sumZ_vec,corrmat3d_gpu,     &
                                        pft1,pft2,                  &
                                        sqsums_refs,sqsums_ptcls,   &
                                        ipart,nrot,nk,              &
                                        lda,ldb,ldc,alpha,          &
                                        s_bench, s_debug_gpu)
     call gettimeofday_c(et_ZM_r)
     call elapsed_time_c(st_ZM_r,et_ZM_r,elps_corr_ZM_gpu)
     call stop_timer_cpu("corr_ZM_gpu")
     if (err /= RC_SUCCESS ) write(*,*) "get_polarft_multi_gpus_gpu_c_ = -1"
     !$omp parallel default(shared) private(iptcl)
     !$omp workshare
     corrmat2d_gpu = maxval(corrmat3d_gpu, dim=3)
     inplmat2d_gpu = maxloc(corrmat3d_gpu, dim=3)
     !$omp end workshare
     !$omp end parallel

     where ( corrmat2d_gpu ==  infinity ) corrmat2d_gpu = 0.0
     where ( corrmat3d_gpu ==  infinity ) corrmat3d_gpu = 0.0
     where ( corrmat3d_gpu == -infinity ) corrmat3d_gpu = 0.0
     
     write(*,*)"Elapsed time elps_corr_ZM_gpu: ",real(elps_corr_ZM_gpu),"(seconds)"
     write(*,*)'                                                               '
     write(*,*)'***************Comparing cormat3d from CPU vs GPU**************'
     write(*,*)'                                                               '

     write(*,'(10x,a,5x,a,10x,a,4x,a,5x,a,1x,a)') &
          "ipart","irot","ik","(corrmat3d)","(corrmat3d_gpu)"
     do i=1,ipart - (ipart-2)
        do irot=1,nrot - (nrot-2)
           do ik=1,nk - (nk-2)
              write(*,*)i, irot, ik, corrmat3d(i,irot,ik), corrmat3d_gpu(i,irot,ik)
           end do
        end do
     end do

     speedup = real(elps_corr_cpu) / real(elps_corr_ZM_gpu)
     sum_diff_crmt2d = sum(corrmat2d-corrmat2d_gpu)
     sum_diff_crmt3d = sum(corrmat3d-corrmat3d_gpu)
     write(*,*) "Differences sum(corrmat3d-corrmat3d_gpu): ",sum_diff_crmt3d
     write(*,*) "Differences sum(corrmat2d-corrmat2d_gpu): ",sum_diff_crmt2d
     write(*,*) "N       corrmat2d(npart,npart): ",size_crmt2d/1.e6,"(Millions)"
     write(*,*) "N  corrmat3d(npart,npart,nrot): ",size_crmt3d/1.e6,"(Millions)"
     write(*,*)'***************************************************************'
     write(*,'(x,a,i2,6x,a,i3)') &
          " Kernel: ",s_polar%ikrnl,"threads/Block: ",s_polar%threadsPerBlock
     write(*,*) " Speed up from CPU to GPU     : ",speedup
     write(*,*)'***************************************************************'

#if defined (BENCH)
     if( ibench_write ) then
        open(4,file='MltGPU_GnCrA_3D_ZM.asc',status='unknown',position='append')
        write(4,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
             real(elps_corr_cpu) / real(elps_corr_ZM_gpu)
        close(4)
     end if
#endif     

!*******************************************************************************
!    Data deallocation
!
!*******************************************************************************

     !Freeing the CPU memory 
     deallocate(sumZ_vec)
     deallocate(sqsums_refs)
     deallocate(sqsums_ptcls)

     deallocate(pft1)
     deallocate(pft2)

     deallocate(corrmat3d)
     deallocate(hadamard_prod)
     deallocate(corrmat2d)
     deallocate(inplmat2d)
     
     !Freeing the cpu variables used for gpu calculation
     deallocate(corrmat3d_gpu)
     deallocate(corrmat2d_gpu)
     deallocate(inplmat2d_gpu)

  end do !end of do ipart = start_npart, npart, istep
  
  !freeing ressources on host for the data structure 
  deallocate(a_devD)
  deallocate(devname)
!*******************************************************************************
!    Environment shutdown
!
!*******************************************************************************
  
  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_deviceQuery_gpu()

contains
!//////////////////// Helper subroutine and functions //////////////////////////
!//
!*******************************************************************************
!    Subroutine to calculate the correlation matrix
!
!*******************************************************************************
  subroutine get_polarft_gencorrall_cpu(corrmat3d, hadamard_prod, &
                                        d1lim, d2lim,             &
                                        s_polar,                  &
                                        sumZ_vec,                 &
                                        pft1,pft2,                &
                                        sqsums_refs,sqsums_ptcls, &
                                        ipart,nrot,nk,            &
                                        lda,ldb,ldc,alpha)
    !$ use omp_lib
    !$ use omp_lib_kinds
    implicit none
    !global variables
    type(polar_corr_calc) :: s_polar
    integer               :: ipart
    integer               :: nrot
    integer               :: nk

    integer               :: lda
    integer               :: ldb,ldc
    real(sp)              :: alpha

    integer, dimension(2) :: d1lim, d2lim        ! dimension arrays
    real(sp)              ::      sumZ_vec(*)
    real(sp)              ::   sqsums_refs(*)    !self%sqsums_refs(self%nrefs)
    real(sp)              ::  sqsums_ptcls(*)    !self%sqsums_ptcls(self%nptcls)
    complex(sp)           ::          pft1(ipart,nrot/2,*)
    complex(sp)           ::          pft2(2*ipart,2*nrot,*)
    real(sp)              ::     corrmat3d(ipart,ipart,*)
    real(sp)              :: hadamard_prod(ipart,nrot/2,*)
    
    !local variables
    integer               :: refsz
    !indexers
    integer               :: jpart,irot,iptcl
    !indexers temporary
    integer  :: j
    
    !start of the execution commnads

    refsz = nrot / 2
    
    d1lim(2) = 0
    d2lim(2) = 0
    do jpart = 1,ipart
       d1lim(1) = jpart
       d1lim(2) = jpart + ipart - 1
       do irot = 1,nrot
          d2lim(1) = irot 
          d2lim(2) = irot + refsz - 1

          !$omp parallel default(shared) private(iptcl)
          !$omp workshare
          hadamard_prod(:ipart,:nrot/2,:nk) = real(pft1(:ipart,:nrot/2,:nk)*&
               conjg(pft2(d1lim(1):d1lim(2),d2lim(1):d2lim(2),:nk)))
          !$omp end workshare
          
          ! correlation normalization

          !$omp do schedule(auto)
          do iptcl=1,ipart
             corrmat3d(iptcl,jpart,irot) = sum(hadamard_prod(iptcl,:nrot/2,:nk))
          end do
          !$omp end do

          !$omp workshare
          corrmat3d(:,jpart,irot) = corrmat3d(:,jpart,irot)/&
               sqrt( sqsums_refs(:ipart) * sqsums_ptcls(d1lim(1):d1lim(2)) )
          !$omp end workshare nowait

          !$omp end parallel

       end do
    end do

    write(*,*)d1lim, d2lim
    
    return
  end subroutine get_polarft_gencorrall_cpu
!*******************************************************************************
!    Routine to print the differences between two, ccorrmat3d and hadamard_prod
!    and the differences from the Subroutine to calculate the correlation matrix
!
! real(sp), allocatable ::     corrmat3d(:,:,:)!cor matrix 3d final output
! real(sp), allocatable :: hadamard_prod(:,:,:)!(self%nptcls,self%refsz,self%nk)
!
!    allocate(    corrmat3d(ipart, ipart, nrot) )
!    allocate(hadamard_prod(ipart, refsz, nk) )
!    if (debug) call get_Hadmr_cormat3d_diff(hadamard_prod, hadamard_prod_sbr, &
!                                            corrmat3d, corrmat3d_sbr,         &
!                                            ipart,nrot,nk)
!*******************************************************************************
  subroutine get_Hadmr_cormat3d_diff(hadamard_prod, hadamard_prod_sbr, &
                                     corrmat3d, corrmat3d_sbr, ipart,nrot,nk,nprint)
    implicit none

    !global variables
    integer               :: ipart
    integer               :: nrot
    integer               :: nk
    integer               :: nprint
    real(sp)              ::     corrmat3d(ipart,ipart,*)
    real(sp)              :: hadamard_prod(ipart,nrot/2,*)

    real(sp)              ::     corrmat3d_sbr(ipart,ipart,*)
    real(sp)              :: hadamard_prod_sbr(ipart,nrot/2,*)
    
    !local variables
    
     write(*,'(10x,a,5x,a,10x,a,5x,a,4x,a)')"ipart","irot","ik","(corrmat3d)","(corrmat3d_sbr)"
     do i=1,ipart - (ipart-nprint)
        do irot=1,nrot - (nrot-nprint)
           do ik=1,nk - (nk-nprint)
              write(*,*)i, irot,ik,corrmat3d(i,irot,ik), corrmat3d_sbr(i,irot,ik)
           end do
        end do
     end do
     write(*,'(10x,a,5x,a,10x,a,3x,a,2x,a)')"ipart","irot","ik","(hadamard_prod)","(hadamard_prod_sbr)"
     do i=1,ipart - (ipart-nprint)
        do irot=1,nrot - (nrot-nprint)
           do ik=1,nk - (nk-nprint)
              write(*,*)i, irot,ik,hadamard_prod(i,irot,ik), hadamard_prod_sbr(i,irot,ik)
           end do
        end do
     end do
     
     write(*,*) "Differences sum(corrmat3d -   corrmat3d_sbr): ", &
          sum(corrmat3d(:ipart,:ipart,:nk)-corrmat3d_sbr(:ipart,:ipart,:nk))
     write(*,*) "Differences sum(hadamard_prod - hadamard_prod_sbr): ", &
          sum(hadamard_prod(:ipart,:nrot/2,:nk)-hadamard_prod_sbr(:ipart,:nrot/2,:nk))

    return
  end subroutine get_Hadmr_cormat3d_diff
!*******************************************************************************
!   memoizer for the sqsums_ref  
!
!*******************************************************************************
  subroutine memoize_sqsum_ref( sqsums_refs, pfts_refs, nref, nrot, nk )
    use simple_math, only: csq
    integer :: iref
    integer, intent(in) :: nref, nrot, nk
    real(sp) :: sqsums_refs(*)       !< memoized square sums for the correlation calculations
    complex(sp) :: pfts_refs(nref, nrot, * )

    do iref = 1,nref
       sqsums_refs(iref) = sum( csq(pfts_refs(iref,:nrot,:nk)) )
    end do

    return
  end subroutine memoize_sqsum_ref
!*******************************************************************************
!   memoizer for the sqsums_ptcl
!
!*******************************************************************************
  subroutine memoize_sqsum_ptcl( sqsums_ptcls, pfts_ptcls, two_ipart, refsz, nk)
    use simple_math, only: csq
    integer, intent(in) :: refsz, nk, two_ipart
    integer :: ptcl_ind
    integer :: iptcl
    real(sp) :: sqsums_ptcls(*)
    complex(sp) :: pfts_ptcls(two_ipart, refsz, * )
    !local variable
    integer :: i
    !start of the execution commands
    iptcl = two_ipart / 2
    !$omp parallel default(shared) private(i)
    !$omp do schedule(auto)
    do i = 1,iptcl  !ptcl_ind, two_ipart - ptcl_ind + 1
       ptcl_ind = ptcl_index(iptcl, i)
       sqsums_ptcls(i) = sum(csq(pfts_ptcls(i:ptcl_ind,:refsz,:nk)))
    end do
    !$omp end do
    !$omp end parallel
    return
  end subroutine memoize_sqsum_ptcl
!*******************************************************************************
!  indexer for the particle indxer
!
!*******************************************************************************
  function ptcl_index( two_ipart, iptcl_in ) result( iptcl_out )
    integer, intent(in)                 :: iptcl_in
    integer :: iptcl_out
    integer :: two_ipart
    iptcl_out = two_ipart + iptcl_in - 1
  end function ptcl_index

end program testing_MultiGPU_gpu
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
