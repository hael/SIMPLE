! ============================================================================
! Name        : testing_corr_gpu.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the corr calc calculation with CUDA
!             :
! ============================================================================
!
program testing_corr_gpu
  use, intrinsic :: iso_c_binding
  use simple_defs
 ! use simple_cuda_defs
  use simple_mpi_defs
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
!  use simple_deviceQuery_gpu
  use simple_math, only: calc_corr, csq

  implicit none
!#if defined (SIMPLE_MPI)
!  include 'mpif.h'
!#endif
type, bind(c) :: polar_corr_calc
   real(c_float) :: r_polar
   real(c_float) :: sumasq_polar
   real(c_float) :: sumbsq_polar
end type polar_corr_calc

#define devptr_t integer*8
!  type(deviceQuery_gpu)         :: devQ
!  type(deviceDetails)             :: devD
  type(polar_corr_calc)         :: s_polar
  integer, parameter            :: npart=1200
  integer, parameter            :: start_npart=1200, istep=125
  integer, parameter            :: nrot = 400
  integer, parameter            :: nk = 59

  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc

  !local variables
  integer                       :: ndev

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
  complex(sp), allocatable      :: pft1(:,:,:)
  complex(sp), allocatable      :: pft2(:,:,:)

  !gpu variables for the calculation of the r, sumasq and sumbsq

  real(sp)                      :: r_gpu, sumasq_gpu, sumbsq_gpu
  real(sp)                      :: dble_alpha,dble_beta

  !timer variables

  double precision              :: elps_corr_FN_gpu
  double precision,dimension(2) :: st_FN_r, et_FN_r

  double precision              :: elps_corr_NN_gpu
  double precision,dimension(2) :: st_NN_r, et_NN_r

  double precision              :: elps_corr_cpu
  double precision,dimension(2) :: st_corr_cpu, et_corr_cpu

  !openMP variables
  integer                       :: chunk
  !MPI variables
!#if defined (SIMPLE_MPI)
!  character(MPI_MAX_PROCESSOR_NAME) :: hostname
!#endif
  !counters

  integer                       :: i,jpart,ik,irot
  integer                       :: i1,i2

  integer :: get_polarft_corr_gpu_c_

  !start of the execution commands
  !start of the greeting message
!  call hello_gpu_magma()
  call timestamp()
  call start_Alltimers_cpu()

  !starting the mpi environment
  call hello_mpi(my_rank,n_proc)
  call simple_mpi_init(mpierr,my_rank,n_proc)
  call MPI_GET_PROCESSOR_NAME(hostname, len, mpierr)
  write(*,*)'Number of proc = ',n_proc,' My rank = ',my_rank, &
       "Length of host = ",len,'Running on = ',hostname

  write(*,*)MASTER,FROM_MASTER,FROM_WORKER
  write(*,*)mpierr
  
  !starting the cuda environment
!  call simple_cuda_init(err)
!  if (err .ne. 0 ) write(*,*) 'cublas init failed'

!  call devQ%new_deviceQuery_gpu(devD)

!  ndev = devQ%get_devCnt()
!  write(*,*)'cudaGetDeviceCount returned: ',ndev
  !writting the number of devices from data structure
!  write(*,*)'cudaGetDeviceCount returned from data structure: ',devD%ndev
!  if ( ndev /= devD%ndev) write(*,*)"data structure and getter do not match"

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
     write(*,'(4x,a,4x,a,5x,a,4x,a,3x,a,3x,a,8x,a)') &
          "inpart","nrot","nk","nrad1","nrad2","klp","khp"
     write(*,'(7(3x,i5))') ipart,nrot,nk,nradial1, nradial2, klp, khp
     write(*,*) "Number of elements in pft(npart,nrot,nk): ",ipart*nrot*nk
     write(*,*) "Number of elements in pft(npart,nrot,nk): ",ipart*nrot*nk/1.e6,"(Millions)"
     write(*,*) 

     lda = ipart

     !allocating the complex matrices
     allocate(pft1(ipart,nrot,nk))
     allocate(pft2(ipart,nrot,nk))

     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft1)
     call getCRandomGaussianDistr_3D(ipart,nrot,nk,pft2)

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
     dble_beta = 0.0

     !proceeding to GPU calculation
     r_gpu = 1.0
     sumasq_gpu = 2.0
     sumbsq_gpu = 3.0
     s_polar%r_polar = 1.342
     s_polar%sumasq_polar = 2.132
     s_polar%sumbsq_polar = 3.123
     write(*,*)'                           N    N                              '
     call start_timer_cpu("corr_NN_gpu")
     call gettimeofday_c(st_NN_r)
!     err = get_polarft_corr_gpu_c_(s_polar,"N","N",            &
!                                   r_gpu,sumasq_gpu,sumbsq_gpu,&
!                                   pft1,pft2,                  &
!                                   ipart,nrot,nk,              &
!                                   lda,ldb,ldc,dble_alpha,dble_beta)
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
!     err = get_polarft_corr_gpu_c_(s_polar,"F","N",            &
!                                   r_gpu,sumasq_gpu,sumbsq_gpu,&
!                                   pft1,pft2,                  &
!                                   ipart,nrot,nk,              &
!                                   lda,ldb,ldc,dble_alpha,dble_beta)
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

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'                                                               '
     
     !Freeing the CPU memory 
     deallocate(pft1)
     deallocate(pft2)

#if defined (BENCH)
     open(1,file='corr_calc_NN_.asc',status='unknown',position='append')
     write(1,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
          real(elps_corr_cpu) / real(elps_corr_NN_gpu)
     close(1)

     open(2,file='corr_calc_FN_.asc',status='unknown',position='append')
     write(2,'(1x,f15.8,2x,f15.8)')ipart*nrot*nk/1.e6, &
          real(elps_corr_cpu) / real(elps_corr_FN_gpu)
     close(2)
#endif     

  end do

  !end  of the MPI computation
  call simple_mpi_finalize(mpierr)
  write(*,*)'Printing the result'

  !shutting down the environment
  ! call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
!  call bye_gpu_magma()

end program testing_corr_gpu
