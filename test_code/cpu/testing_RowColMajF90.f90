! ============================================================================
! Name        : testing_RowColMajF90
! Author      : Frederic Bonnet
! Version     :
! Date        : 11th of November 2015
! Description : tests the summing and corr calc calculation
!             :
! ============================================================================
!
program testing_RowColMajF90
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
  integer, parameter            :: npart=3
  integer, parameter            :: start_npart=3, istep=128
  integer, parameter            :: nrot = 3
  integer, parameter            :: nk = 3

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

     deallocate(r1_vec)
     deallocate(r2_vec)

     write(*,*)'                                                               '
     write(*,*)'***************indexer checker Row and Colum major*************'
     write(*,*)'                                                               '

     hadamard_prod = 0.0
     call get1to9ColMajDspMat_3D_cpu(ipart,nrot,nk,lda,hadamard_prod)

     do jpart=1,ipart
        do irot=1,nrot
           do ik=1,nk
              write(*,*) jpart,irot,ik,hadamard_prod(jpart,irot,ik)
           end do
        end do
     end do



     !Freeing the CPU memory 
     deallocate(pft1)
     deallocate(pft2)
     deallocate(hadamard_prod)

  end do

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_RowColMajF90
