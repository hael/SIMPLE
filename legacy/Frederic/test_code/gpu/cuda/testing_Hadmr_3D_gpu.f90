! ============================================================================
! Name        : testing_Hadmr_3D_gpu
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests various kernel for Hadamard product and summing in CUDA
!             :
! ============================================================================
!
program testing_Hadmr_3D_gpu
  use, intrinsic :: iso_c_binding
  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_deviceQuery_gpu
  use simple_math, only: calc_corr, csq
  implicit none
  logical, parameter :: ibench_local=.true.       !< benchmark result or not
  logical, parameter :: ibench_write_local=.true. !write benchmark result or not

  logical, parameter :: debug=.false.         !< debug indicator
  logical, parameter :: debug_cpu=.false.     !false/true
  logical, parameter :: debug_high=.true.     !true/false
  logical, parameter :: debug_write=.false.   !false/true
  logical, parameter :: debug_write_C=.false. !false/true

  !openMP variables
  integer                       :: nthr=8

#define devptr_t integer*8
  type(deviceQuery_gpu)         :: devQ
  type(deviceDetails)           :: devD
  type(polar_corr_calc)         :: s_polar
  type(t_debug_gpu)             :: s_debug_gpu
  type(t_bench)                 :: s_bench
  integer, parameter            :: npart=16
  integer, parameter            :: start_npart=16, istep=1
  integer, parameter            :: nrot = 200
  integer, parameter            :: nk = 100

  integer                       :: ikernel = 6           !Kernel option
  integer                       :: threadsPerBlock = 512 !threadsPerBlock
  integer                       :: nx_3D=16, ny_3D=16, nz_3D=4 !3D mesh
  
  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc

  !local variables
  type(deviceDetails),allocatable :: a_devD(:)
  !gpu gear
  integer                       :: ndev
  integer                       :: rc

  !index variable major
  integer                       :: ipart

  !variables for the function 

  integer                       :: nradial1    !< #radial vectors (angle index) 
  integer                       :: nradial2    !< #radial vectors (angle index)
  integer                       :: klp         !< low-pass frequency limit
  integer                       :: khp         !< high-pass frequency limit 
  integer                       :: rot1, rot2

  real(sp)                      :: alpha,beta  
  real(sp)                      :: r, sumasq, sumbsq
  real(sp), allocatable         :: r1_vec(:)
  real(sp), allocatable         :: r2_vec(:)
  real(sp), allocatable         :: sumP_vec(:)
  real(sp), allocatable         :: sumX_vec(:)
  real(sp), allocatable         :: hadamard_prod(:,:,:)
  real(sp), allocatable         :: hadamard_prod_1D(:)
  real(sp), allocatable         :: hadamard_prod_23D(:,:)
  complex(sp), allocatable      :: pft1(:,:,:)
  complex(sp), allocatable      :: pft2(:,:,:)

  !gpu variables for the calculation of the r, sumasq and sumbsq

  real(sp)                      :: r_gpu, sumasq_gpu, sumbsq_gpu
  real(sp)                      :: dble_alpha,dble_beta

  !timer variables

  double precision              :: elps_corr_XN_gpu
  double precision,dimension(2) :: st_XN_r, et_XN_r

  double precision              :: elps_corr_PN_gpu
  double precision,dimension(2) :: st_PN_r, et_PN_r

  double precision              :: elps_corr_FN_gpu
  double precision,dimension(2) :: st_FN_r, et_FN_r

  double precision              :: elps_corr_NN_gpu
  double precision,dimension(2) :: st_NN_r, et_NN_r

  double precision              :: elps_corr_cpu
  double precision,dimension(2) :: st_corr_cpu, et_corr_cpu

  !indexers
  integer                       :: i,jpart,ik,irot
  integer                       :: i1,i2
  integer                       :: idev

  !functions calls
  integer :: get_dev_count_c
  integer :: get_polarft_corr_gpu_c

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
  !start of the greeting message
  call hello_gpu_magma()
  call timestamp()
  call start_Alltimers_cpu()

!$ call omp_set_num_threads(nthr)

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. 0 ) write(*,*) 'cublas init failed'

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

  write(*,*)'cudaGetDeviceCount returned: ',ndev

  do ipart = start_npart, npart, istep

     lda = ipart
     ldb = nrot
     ldc = ipart
!*******************************************************************************
!     now testign the corr function 
!
!*******************************************************************************

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     now testign the corr function                             '
     write(*,*)'***************************************************************'

     rot1 = 1
     rot2 = 1
     nradial1 = ipart
     nradial2 = ipart
     klp = ipart
     khp = 1

     write(*,*) 
     write(*,'(4x,a,3x,a,5x,a,4x,a,3x,a,3x,a,8x,a)') &
          "inpart","nrot","nk","nrad1","nrad2","klp","khp"
     write(*,'(7(3x,i5))') ipart,nrot,nk,nradial1, nradial2, klp, khp
     write(*,*) "Number of elements in pft(npart,nrot,nk): ",ipart*nrot*nk
     write(*,*) "Number of elements in pft(npart,nrot,nk): ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*) 

     !allocating the complex matrices
     allocate(pft1(ipart,nrot,nk))
     allocate(pft2(ipart,nrot,nk))
     allocate(sumP_vec(ipart))
     allocate(sumX_vec(ipart))

     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft1)
     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft2)
!     call getCSeq_3D(ipart,nrot,nk,pft1)
!     call getCSeq_3D(ipart,nrot,nk,pft2)

     write(*,'(10x,a,5x,a,10x,a,16x,a,26x,a)')"ipart","irot","ik","(pft1)","(pft2)"
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
                   real( pft1(i,irot,ik)*conjg(pft2(i,irot,ik)) ), &
                   real(pft1(i,irot,ik)) * real(pft2(i,irot,ik)) + &
                   imag(pft1(i,irot,ik)) * imag(pft2(i,irot,ik)) 
           end do
        end do
     end do

     write(*,*)'                                                               '
     write(*,*)'***************CPU corr****************************************'
     write(*,*)'                                                               '

     i1     = rot1
     i2     = rot2
     sumasq = 0.
     sumbsq = 0.
     r      = 0.
     call start_timer_cpu("corr_cpu")
     !  do i=1,nradial1/2

     call gettimeofday_c(st_corr_cpu)

     open(1,file='AAstar_cpu_GNUf.log',status='unknown',action='write')

     do jpart=1,ipart
        do irot=1,nrot
           do ik=1,nk
              r = r+real(pft1(i1,irot,ik)*conjg(pft2(i2,irot,ik)))
              sumasq = sumasq+csq(pft1(i1,irot,ik))
              sumbsq = sumbsq+csq(pft2(i2,irot,ik))
!              write(1,'(i5,x,i5,x,i5,x,f20.8)')i1, irot, ik,csq(pft1(i1,irot,ik))
           end do
        end do
        i1 = i1+1
        if( i1 > nradial1 ) i1 = 1
        i2 = i2+1
        if( i2 > nradial2 ) i2 = 1
        !write(*,*)i1,i2,i,j,r
     end do

     close(1)

     call gettimeofday_c(et_corr_cpu)
     call elapsed_time_c(st_corr_cpu,et_corr_cpu,elps_corr_cpu)

     call stop_timer_cpu("corr_cpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     !write(*,*)"Full sum of real(pft1) is: ",sum(real(pft1))
     write(*,*)"Full sum of r is         : ",r
     write(*,*)"Full sum of sumasq is    : ",sumasq
     write(*,*)"Full sum of sumbsq is    : ",sumbsq
     write(*,*)"the correlator"
     r = calc_corr(r,sumasq*sumbsq)
     write(*,*)"after calc_corr r        : ",r
     write(*,*)"Elapsed time for corr_cpu: ",real(elps_corr_cpu),"(seconds)"

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     dble_alpha = 1.0
     r_gpu = 1.0
     sumasq_gpu = 2.0
     sumbsq_gpu = 3.0
     s_polar%r_polar = 1.342
     s_polar%sumasq_polar = 2.132
     s_polar%sumbsq_polar = 3.123
     s_polar%ikrnl           = ikernel
     s_polar%threadsPerBlock = threadsPerBlock
     s_polar%nx              = nx_3D
     s_polar%ny              = ny_3D
     s_polar%nz              = nz_3D     
     write(*,*)'                           N    N                              '
     call start_timer_cpu("corr_NN_gpu")
     call gettimeofday_c(st_NN_r)
     err = get_polarft_corr_gpu_c(a_devD,                     &
                                  s_polar,"N","N",            &
                                  r_gpu,                      &
                                  pft1,pft2,                  &
                                  ipart,nrot,nk,              &
                                  lda,ldb,ldc,dble_alpha,     &
                                  s_bench, s_debug_gpu)
     call gettimeofday_c(et_NN_r)
     call elapsed_time_c(st_NN_r,et_NN_r,elps_corr_NN_gpu)
     call stop_timer_cpu("corr_NN_gpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Full sum of r is         : ",s_polar%r_polar
     write(*,*)"Full sum of sumasq is    : ",s_polar%sumasq_polar
     write(*,*)"Full sum of sumbsq is    : ",s_polar%sumbsq_polar
     write(*,*)"the correlator"
     r_gpu = calc_corr(s_polar%r_polar,s_polar%sumasq_polar*s_polar%sumbsq_polar)
     write(*,*)"after calc_corr r_gpu    : ",r_gpu
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_NN_gpu),"(seconds)"

     write(*,*)'                           F    N                              '

     call start_timer_cpu("corr_FN_gpu")
     call gettimeofday_c(st_FN_r)
     err = get_polarft_corr_gpu_c(a_devD,                     &
                                  s_polar,"F","N",            &
                                  r_gpu,                      &
                                  pft1,pft2,                  &
                                  ipart,nrot,nk,              &
                                  lda,ldb,ldc,dble_alpha,     &
                                  s_bench, s_debug_gpu)
     call gettimeofday_c(et_FN_r)
     call elapsed_time_c(st_FN_r,et_FN_r,elps_corr_FN_gpu)
     call stop_timer_cpu("corr_FN_gpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Full sum of r is         : ",s_polar%r_polar
     write(*,*)"Full sum of sumasq is    : ",s_polar%sumasq_polar
     write(*,*)"Full sum of sumbsq is    : ",s_polar%sumbsq_polar
     write(*,*)"the correlator"
     r_gpu = calc_corr(s_polar%r_polar,s_polar%sumasq_polar*s_polar%sumbsq_polar)
     write(*,*)"after calc_corr r_gpu    : ",r_gpu
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_FN_gpu),"(seconds)"

     write(*,*)'                           P    N                              '

     call start_timer_cpu("corrPN_gpu")
     call gettimeofday_c(st_PN_r)
     err = get_polarft_corr_gpu_c(a_devD,                     &
                                  s_polar,"P","N",            &
                                  sumP_vec,                   &
                                  pft1,pft2,                  &
                                  ipart,nrot,nk,              &
                                  lda,ldb,ldc,dble_alpha,     &
                                  s_bench, s_debug_gpu)
     call gettimeofday_c(et_PN_r)
     call elapsed_time_c(st_PN_r,et_PN_r,elps_corr_PN_gpu)
     call stop_timer_cpu("corrPN_gpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_PN_gpu),"(seconds)"

     write(*,*)'                           X    N                              '

     call start_timer_cpu("XN_gpu")
     call gettimeofday_c(st_XN_r)
     err = get_polarft_corr_gpu_c(a_devD,                     &
                                  s_polar,"X","N",            &
                                  sumX_vec,                   &
                                  pft1,pft2,                  &
                                  ipart,nrot,nk,              &
                                  lda,ldb,ldc,dble_alpha,     &
                                  s_bench, s_debug_gpu)
     call gettimeofday_c(et_XN_r)
     call elapsed_time_c(st_XN_r,et_XN_r,elps_corr_XN_gpu)
     call stop_timer_cpu("XN_gpu")

     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*)"Elapsed time for corr_gpu: ",real(elps_corr_XN_gpu),"(seconds)"

     write(*,*)'                                                               '
     write(*,*)'***************CPU corr double sum 2 and 3D********************'
     write(*,*)'                                                               '

     !allocating the complex matrices
     allocate(r1_vec(ipart))
     allocate(r2_vec(ipart))
     allocate(hadamard_prod(ipart,nrot,nk))

     allocate(hadamard_prod_1D(ipart))
     allocate(hadamard_prod_23D(nrot,nk))

     hadamard_prod(:,:,:) = real(pft1(:,:,:)*conjg(pft2(:,:,:)))
     r1_vec = 0.0
     r2_vec = 0.0
     do jpart=1,ipart
        r1_vec(jpart) = sum(hadamard_prod(jpart,:,:))
        do irot=1,nrot
           do ik=1,nk
              r2_vec(jpart) = r2_vec(jpart) + hadamard_prod(jpart,irot,ik)
              hadamard_prod_1D(jpart) = hadamard_prod(jpart,irot,ik)
              hadamard_prod_23D(irot,ik) = hadamard_prod(jpart,irot,ik)
           end do
        end do
     end do

     write(*,'(11x,a,3x,a,11x,a,11x,a,6x,a)')"i","r1_vec","r2_vec","sumP_vec[i]","sumX_vec[i]"
     do i=1,ipart - (ipart-5)
        write(*,*)i, r1_vec(i), r2_vec(i),sumP_vec(i),sumX_vec(i)
     end do
     write(*,*) "Differences sum(r1_vec -   r2_vec): ",sum(r1_vec-r2_vec)
     write(*,*) "Differences sum(r1_vec - sumP_vec): ",sum(r1_vec-sumP_vec)
     write(*,*) "Differences sum(r1_vec - sumX_vec): ",sum(r1_vec-sumX_vec)
     write(*,*) "Relative Differences sum((r1_vec -   r2_vec)/r1_vec): ",sum((r1_vec-r2_vec)/r1_vec)
     write(*,*) "Relative Differences sum((r1_vec - sumP_vec)/r1_vec): ",sum((r1_vec-sumP_vec)/r1_vec)
     write(*,*) "Relative Differences sum((r1_vec - sumX_vec)/r1_vec): ",sum((r1_vec-sumX_vec)/r1_vec)
     write(*,*)"size(pft(npart,nrot,nk)) : ",ipart*nrot*nk/1.e6,"(Millions)"

     deallocate(hadamard_prod)
     deallocate(hadamard_prod_1D)
     deallocate(hadamard_prod_23D)

     deallocate(r1_vec)
     deallocate(r2_vec)
     deallocate(sumP_vec)
     deallocate(sumX_vec)

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'

     !Freeing the CPU memory 
     deallocate(pft1)
     deallocate(pft2)

#if defined (BENCH)
     if( ibench_write ) then
        open(1,file='corr_calc_NN.asc',status='unknown',position='append')
        write(1,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
             real(elps_corr_cpu) / real(elps_corr_NN_gpu)
        close(1)

        open(2,file='corr_calc_FN.asc',status='unknown',position='append')
        write(2,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
             real(elps_corr_cpu) / real(elps_corr_FN_gpu)
        close(2)
        
        open(3,file='corr_calc_PN.asc',status='unknown',position='append')
        write(3,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
             real(elps_corr_cpu) / real(elps_corr_PN_gpu)
        close(3)
     
        open(4,file='corr_calc_XN.asc',status='unknown',position='append')
        write(4,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
             real(elps_corr_cpu) / real(elps_corr_XN_gpu)
        close(4)
     end if
#endif     

  end do

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_Hadmr_3D_gpu
