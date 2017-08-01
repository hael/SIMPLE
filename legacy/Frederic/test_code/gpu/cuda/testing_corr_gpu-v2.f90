! ============================================================================
! Name        : timer_tester.f90
! Author      : Frederic Bonnet
! Version     :
! Date        : 30th of March 2015
! Description : tests the functionality of the My_zgetri on cuda
!             :
! ============================================================================
!
program testing_corr_gpu

  use simple_defs
  use simple_cuda_defs
  use simple_cuda
  use matrixGetter
  use invert_gpu_utils
  use greeting_version
  use simple_timing
  use simple_math, only: calc_corr_dble, csq_dble

  implicit none
#define devptr_t integer*8
  integer, parameter            :: nmx=10541
  integer, parameter            :: start_nmx=10541, istep=1000

  integer                       :: err
  integer                       :: size_of_elt
  integer                       :: lda, n
  integer                       :: ldb,ldc
  integer                       :: k

  !timing variables

  integer                       :: inmx

  !index variable
  integer                       :: ilda
  integer                       :: jlda

  !variables for the function 

  integer                       :: nradial1    !< #radial vectors (angle index) (= # of components in each shell)
  integer                       :: nradial2    !< #radial vectors (angle index) (= # of components in each shell)
  integer                       :: klp         !< low-pass frequency limit
  integer                       :: khp         !< high-pass frequency limit 

  integer                       :: rot1, rot2
  integer                       :: ring2

  real(dp)                      :: r, sumasq, sumbsq
  complex(dp), allocatable      :: pft1(:,:)
  complex(dp), allocatable      :: pft2(:,:)

  !gpu variables for the calculation of the r, sumasq and sumbsq

  real(dp)                      :: r_gpu, sumasq_gpu, sumbsq_gpu
  real(dp)                      :: dble_alpha,dble_beta
  real(dp), allocatable         :: DA_gpu(:,:)

  !the pointers devices
  devptr_t                      ::  devPtrA_D1
  devptr_t                      ::  devPtrA_D2
  devptr_t                      ::  devPtrA_D3
  devptr_t                      ::  devPtrA_pft1
  devptr_t                      ::  devPtrA_pft2

  !counters

  integer                       :: i,j,in
  integer                       :: i1,i2

  !timer variables

  double precision              :: elps_sumasq
  double precision,dimension(2) :: st_sumasq, et_sumasq

  double precision              :: elps_sumbsq
  double precision,dimension(2) :: st_sumbsq, et_sumbsq

  double precision              :: elps_r
  double precision,dimension(2) :: st_r, et_r
  double precision              :: elps_sum
  double precision,dimension(2) :: st_sum, et_sum

  !start of the execution commands
  !start of the greeting message
  call hello_gpu_magma()
  call timestamp()
  call start_Alltimers_cpu()

  !starting the cuda environment
  call simple_cuda_init(err)
  if (err .ne. 0 ) write(*,*) 'cublas init failed'

  do inmx = start_nmx, nmx, istep
     print *, inmx

     lda = inmx

!*******************************************************************************
!     now testign the corr function 
!
!*******************************************************************************

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'     now testign the corr function                             '
     write(*,*)'***************************************************************'

     ring2 = 2
 
     rot1 = 1
     rot2 = 1
     nradial1 = inmx
     nradial2 = inmx
     klp = inmx
     khp = 1

     write(*,*) 
     write(*,'(7(5x,a))') "ring2","rot1","rot2","nradial1","nradial2","klp","khp"
     write(*,'(7(7x,i5))') ring2, rot1, rot2, nradial1, nradial2, klp, khp
     write(*,*) 

     lda = inmx

     !allocating the complex matrices
     allocate(pft1(inmx,inmx))
     allocate(pft2(inmx,inmx))

     call get1to9RowMajZdpMat_cpu(inmx,lda,pft1)
     call get1to9ColMajZdpMat_cpu(inmx,lda,pft2)
     
     write(*,'(34x,a,32x,a)')"(pft1) row major","(pft2) column major"
     do i=1,inmx - (inmx-2)
        do j=1,inmx - (inmx-2)
           write(*,*)i, j,pft1(i,j), pft2(i,j)
        end do
     end do
     write(*,'(34x,a)')"conjg(pft2)"
     do i=1,inmx - (inmx-2)
        do j=1,inmx - (inmx-2)
           write(*,*)i, j,conjg(pft2(i,j))
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
     write(*,*)"the correlator"
     !  do i=1,nradial1/2
     do i=1,nradial1
        do j=khp,klp
           r = r+real(pft1(i1,j)*conjg(pft2(i2,j)))
           sumasq = sumasq+csq_dble(pft1(i1,j))
           sumbsq = sumbsq+csq_dble(pft2(i2,j))
           !write(*,*)i, j, r
        end do
        i1 = i1+1
        if( i1 > nradial1 ) i1 = 1
        i2 = i2+1
        if( i2 > nradial2 ) i2 = 1
        !write(*,*)i1,i2,i,j,r
     end do
     call stop_timer_cpu("corr_cpu")

     !write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double  precision matrix r is: ",r
     r = calc_corr_dble(r,sumasq*sumbsq)
     write(*,*)"after calc_corr_dble r= ",r!," sumasq= ",sumasq," sumbsq= ",sumbsq

     write(*,*)'                                                               '
     write(*,*)'***************GPU corr****************************************'
     write(*,*)'                                                               '

     !proceeding to GPU calculation
     dble_alpha = 1.0d0
     dble_beta = 0.0d0
     ldb = lda
     ldc = lda
     !allocating memory
     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft1)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft2)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     !setting up the pft 1 and 2 matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft1, inmx, devPtrA_pft1, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft2, inmx, devPtrA_pft2, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     !Now calluating the r value
     allocate(DA_gpu(inmx,inmx))
     write(*,*) 'before the D1 alloc'
     err = cublas_alloc(inmx*inmx, size_of_double, devPtrA_D1)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     write(*,*) 'after the D1 alloc'
     write(*,*) 'before the D1 set'
     err = cublas_set_matrix (inmx, inmx, size_of_double, DA_gpu, inmx, devPtrA_D1, inmx )
     write(*,*) 'after the D1 set'
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     call start_timer_cpu("r_gpu")
     call gettimeofday_c(st_r)
     call zz2dgemm_ElmtWs_tesla_gpu("N", "N", &
                                    inmx, inmx, inmx,  &
                                    dble_alpha,  &
                                    devPtrA_pft1, lda,  &
                                    devPtrA_pft2, ldb,  &
                                    dble_beta,  &
                                    devPtrA_D1, ldc)
     call gettimeofday_c(et_r)
     call elapsed_time_c(st_r,et_r,elps_r)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_r(s)= ",elps_r
     call stop_timer_cpu("r_gpu")

     DA_gpu = 0.0
     write(*,*) 'before the D1 get'
     err = cublas_get_matrix ( inmx, inmx, size_of_double, devPtrA_D1, inmx, DA_gpu, inmx)
     write(*,*) 'after the D1 get'
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     r_gpu = sum(DA_gpu)
     !freein the ressources on device for the first r calculation.
     err = cublas_free(devPtrA_D1)
     deallocate(DA_gpu)

     !***********Now calculate the sumsq(a and b)*********************
     !******sumasq*****
     allocate(DA_gpu(inmx,inmx))

     write(*,*) 'before the D2 alloc'
     err = cublas_alloc(inmx*inmx, size_of_double, devPtrA_D2)
     write(*,*) 'after the D2 alloc'
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double, DA_gpu, inmx, devPtrA_D2, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     call gettimeofday_c(st_sumasq)
     call zz2dgemm_ElmtWs_tesla_sumsq_gpu("N", "N", &
                                          inmx, inmx, inmx,  &
                                          dble_alpha,  &
                                          devPtrA_pft1, lda,  &
                                          devPtrA_pft2, ldb,  &
                                          dble_beta,  &
                                          devPtrA_D2, ldc)
     call gettimeofday_c(et_sumasq)
     call elapsed_time_c(st_sumasq,et_sumasq,elps_sumasq)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_sumasq(s)= ",elps_sumasq

     DA_gpu = 0.0
     write(*,*) 'before the D2 get'
     err = cublas_get_matrix ( inmx, inmx, size_of_double, devPtrA_D2, inmx, DA_gpu, inmx)
     write(*,*) 'after the D2 get'
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     sumasq_gpu = 0.0d0
     sumasq_gpu = sum(DA_gpu)
     err = cublas_free(devPtrA_D2)
     deallocate(DA_gpu)

     !******sumbsq*****
     allocate(DA_gpu(inmx,inmx))
     err = cublas_alloc(inmx*inmx, size_of_double, devPtrA_D3)
     err = cublas_set_matrix (inmx, inmx, size_of_double, DA_gpu, inmx, devPtrA_D3, inmx )

     call gettimeofday_c(st_sumbsq)
     call zz2dgemm_ElmtWs_tesla_sumsq_gpu("N", "N", &
                                          inmx, inmx, inmx,  &
                                          dble_alpha,  &
                                          devPtrA_pft1, lda,  &
                                          devPtrA_pft2, ldb,  &
                                          dble_beta,  &
                                          devPtrA_D3, ldc)
     call gettimeofday_c(et_sumbsq)
     call elapsed_time_c(st_sumbsq,et_sumbsq,elps_sumbsq)
     write(*,'(x,a,x,f20.8)') "the elapse time for the tesla_sumbsq(s)= ",elps_sumbsq

     DA_gpu = 0.0
     err = cublas_get_matrix ( inmx, inmx, size_of_double, devPtrA_D3, inmx, DA_gpu, inmx)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     sumbsq_gpu = 0.0d0
     sumbsq_gpu = sum(DA_gpu)
     err = cublas_free(devPtrA_D3)
     deallocate(DA_gpu)

     !******Now calculating the correlator*****
     r = 0.0d0
     r = calc_corr_dble(r_gpu,sumasq_gpu*sumbsq_gpu)
     write(*,*)"after calc_corr_dble r= ",r!," sumasq= ",sumasq," sumbsq= ",sumbsq

     write(*,*)'                                                               '
     write(*,*)'***************************************************************'
     write(*,*)'                                                               '

     !freeing the ressources on GPU
     err = cublas_free(devPtrA_pft1)
     err = cublas_free(devPtrA_pft2)

     !Freeing the CPU memory 
     deallocate(pft1)  
     deallocate(pft2)

  end do

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_corr_gpu
