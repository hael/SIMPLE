! ============================================================================
! Name        : testing_summing.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 11th of November 2015
! Description : tests the summing and corr calc calculation
!             :
! ============================================================================
!
program testing_summing
  use, intrinsic :: iso_c_binding
  use simple_defs
  use matrixGetter
  use simple_testfunction
  use greeting_version
  use simple_timing
  use simple_systemQuery_cpu
  use simple_math, only: calc_corr, csq

  implicit none

#define devptr_t integer*8
  type(systemQuery_cpu)         :: sysQ
  type(systemDetails)           :: hstD
  integer, parameter            :: npart=200
  integer, parameter            :: start_npart=200, istep=128
  integer, parameter            :: nrot = 40
  integer, parameter            :: nk = 59

  integer                       :: err
  integer                       :: lda
  integer                       :: ldb,ldc

  !local variables
  integer                       :: nCPU_cores
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
  real(sp), allocatable         :: r1_vec(:)
  real(sp), allocatable         :: r2_vec(:)
  real(sp), allocatable         :: hadamard_prod(:,:,:)
  complex(sp), allocatable      :: pft1(:,:,:)
  complex(sp), allocatable      :: pft2(:,:,:)

  !gpu variables for the calculation of the r, sumasq and sumbsq

  real(sp)                      :: r_gpu, sumasq_gpu, sumbsq_gpu
  real(sp)                      :: dble_alpha,dble_beta

  !timer variables

  double precision              :: elps_corr_PN_gpu
  double precision,dimension(2) :: st_PN_r, et_PN_r

  double precision              :: elps_corr_FN_gpu
  double precision,dimension(2) :: st_FN_r, et_FN_r

  double precision              :: elps_corr_NN_gpu
  double precision,dimension(2) :: st_NN_r, et_NN_r

  double precision              :: elps_corr_cpu
  double precision,dimension(2) :: st_corr_cpu, et_corr_cpu

  !openMP variables
  integer                       :: chunk
  !counters

  integer                       :: i,jpart,ik,irot
  integer                       :: i1,i2

  !start of the execution commands
  !start of the greeting message
  call hello_gpu_magma()
  call timestamp()
  call start_Alltimers_cpu()

  call sysQ%new_systemQuery_cpu(hstD)

  write(*,*)'                                                                  '
  write(*,*)'******************************************************************'
  write(*,*)'   Sanity checks on the object(sysQ) and the data structure(hstD) '
  write(*,*)'******************************************************************'
  nCPU_cores = sysQ%get_ncpu_cores()
  write(*,*)'Number of cores on the system                    : ',nCPU_cores
  !writting the number of cores avaible on the system
  write(*,*)'Number of cores returned from data structure hstD: ',hstD%nCPUcores
  if ( nCPU_cores /= hstD%nCPUcores) call sysQ%get_warning_dataStruct()
  write(*,*)'******************************************************************'

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
     write(*,*)'***************CPU corr double summ 2 and 3D*******************'
     write(*,*)'                                                               '

     !allocating the complex matrices
     allocate(r1_vec(ipart))
     allocate(r2_vec(ipart))
     allocate(hadamard_prod(ipart,nrot,nk))

     hadamard_prod(:,:,:) = real(pft1(:,:,:)*conjg(pft2(:,:,:)))
     r1_vec = 0.0
     r2_vec = 0.0
     do jpart=1,ipart
        r1_vec(jpart) = sum(hadamard_prod(jpart,:,:))
        do irot=1,nrot
           do ik=1,nk
              r2_vec(jpart) = r2_vec(jpart) + hadamard_prod(jpart,irot,ik)
           end do
        end do
     end do

     write(*,'(x,a,12x,a,12x,a)')"i","r1_vec","r2_vec"
     do i=1,ipart - (ipart-5)
        write(*,*)i, r1_vec(i), r2_vec(i)
     end do
     write(*,*) "the difference between sum(r1_vec-r2_vec): ",sum(r1_vec-r2_vec)

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'                                                               '

     !Freeing the CPU memory 
     deallocate(pft1)
     deallocate(pft2)
     deallocate(hadamard_prod)
     deallocate(r1_vec)
     deallocate(r2_vec)

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

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_summing
