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
  integer, parameter              :: nmx=1541
  integer, parameter              :: start_nmx=1541, istep=1000

  integer                         :: err
  integer                         :: size_of_elt
  integer                         :: lda, n
  integer                         :: ldb,ldc
  integer                         :: k

  devptr_t                        :: devPtrA_new
  !checking the matrix inversion has work correctly
  devptr_t                        :: devPtrA_invSymZA
  devptr_t                        :: devPtrA_SymZE_gpu
  devptr_t                        :: devPtrA_SymZF_gpu
  !pointers for the correlator kernel
  devptr_t                        ::  devPtrA_sumasq
  devptr_t                        ::  devPtrA_sumbsq
  devptr_t                        ::  devPtrA_DA
  devptr_t                        ::  devPtrA_pft1
  devptr_t                        ::  devPtrA_pft2

  complex(dp)                     :: alpha,beta
  real(dp)                        :: dble_alpha,dble_beta

  real(dp),allocatable            :: matSymZE_gpu(:,:)
  real(dp),allocatable            :: invSymZA(:,:)  !inverted sym Z matrix

  real(dp),allocatable            :: matSymZF_gpu(:,:)  !matrix multiplication 

  !timing variables

  integer                         :: inmx

  !index variable
  integer                         :: ilda
  integer                         :: jlda

  !variables for the function 

  integer                      :: nradial1           !< #radial vectors (angle index) (= # of components in each shell)
  integer                      :: nradial2           !< #radial vectors (angle index) (= # of components in each shell)
  integer                      :: klp               !< low-pass frequency limit
  integer                      :: khp               !< high-pass frequency limit 

  integer                      :: rot1, rot2
  integer                      :: ring2
  !  integer, parameter           :: n=3

  real(dp)                     :: r, sumasq, sumbsq
  real(dp)                     :: r_gpu, sumasq_gpu, sumbsq_gpu
  complex(dp), allocatable     :: pft1(:,:)
  complex(dp), allocatable     :: pft2(:,:)

  real(dp), allocatable        :: DA_gpu(:,:)
  complex(dp), allocatable     :: pft_gpu(:,:)

  double precision,allocatable :: re_matZA(:,:)
  double precision,allocatable :: im_matZA(:,:)

  !counters

  real(dp)                     :: sumDZZ
  integer                      :: i,j,in
  integer                      :: i1,i2

  !timer variables

  double precision                :: elps_corr
  double precision,dimension(2)   :: st_corr, et_corr

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
     write(*,'(7(7x,i3))') ring2, rot1, rot2, nradial1, nradial2, klp, khp
     write(*,*) 

     lda = inmx

     !allocating the complex matrices
     allocate(pft1(inmx,inmx))
     allocate(pft2(inmx,inmx))
     allocate(pft_gpu(inmx,inmx))
     allocate(DA_gpu(inmx,inmx))

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

     call getSymDRandomMat_cpu(inmx,lda,DA_gpu)

     write(*,*)
     write(*,*) "The  DA_gpu"
     write(*,*)
     do ilda = 1,LDA - (inmx-2)
        do jlda = 1,LDA - (inmx-2)
           write(*,*) ilda, jlda, DA_gpu(ilda,jlda)
        end do
     end do
     write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double Complex precision matrix  Real(invSymZA) is: ",sum(DA_gpu)

     write(*,*)'                                                               '
     write(*,*)'***************CPU corr****************************************'
     write(*,*)'                                                               '

     i1     = rot1
     i2     = rot2
     sumasq = 0.
     sumbsq = 0.
     r      = 0.
     call gettimeofday_c(st_corr)
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
     call gettimeofday_c(et_corr)
     call elapsed_time_c(st_corr,et_corr,elps_corr)

     write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double  precision matrix r is: ",r
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
     err = cublas_alloc(inmx*inmx, size_of_double_complex, devPtrA_pft2)
     err = cublas_alloc(inmx*inmx, size_of_double, devPtrA_DA)
     err = cublas_alloc(1, size_of_double, devPtrA_sumasq)
     err = cublas_alloc(1, size_of_double, devPtrA_sumbsq)
     !setting up the matrix on device
     err = cublas_set_matrix (inmx, inmx, size_of_double, DA_gpu, inmx, devPtrA_DA, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft1, inmx, devPtrA_pft1, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (inmx, inmx, size_of_double_complex, pft2, inmx, devPtrA_pft2, inmx )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     sumasq_gpu = 0.0d0
     sumbsq_gpu = 0.0d0
     err = cublas_set_matrix (1, 1, size_of_double, sumasq_gpu, 1, devPtrA_sumasq, 1 )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_set_matrix (1, 1, size_of_double, sumbsq_gpu, 1, devPtrA_sumbsq, 1 )
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     call zz2dgemm_ElmtWs_tesla_gpu("N", "N", &
                                    inmx, inmx, inmx,  &
                                    dble_alpha,  &
                                    devPtrA_pft1, lda,  &
                                    devPtrA_pft2, ldb,  &
                                    dble_beta,  &
                                    devPtrA_DA, ldc, &
                                    devPtrA_sumasq, devPtrA_sumbsq)

     pft_gpu = 0.0
     err = cublas_get_matrix ( inmx, inmx, size_of_double_complex, devPtrA_pft2, inmx, pft_gpu, inmx)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     !write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double Complex precision matrix pft_gpu is: ",sum(pft_gpu)
     DA_gpu = 0.0
     err = cublas_get_matrix ( inmx, inmx, size_of_double, devPtrA_DA, inmx, DA_gpu, inmx)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_get_matrix ( 1, 1, size_of_double, devPtrA_sumasq, 1, sumasq_gpu, 1)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
     err = cublas_get_matrix ( 1, 1, size_of_double, devPtrA_sumbsq, 1, sumbsq_gpu, 1)
     if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

     write(*,*)
     write(*,*) "The DA_gpu(ilda,jlda),"
     write(*,*)
     do ilda = 1,LDA - (inmx-2)
        do jlda = 1,LDA - (inmx-2)
           write(*,*) ilda, jlda, DA_gpu(ilda,jlda), &
                real( pft1(ilda,jlda)*conjg(pft2(ilda,jlda) ) )
        end do
     end do

!     sumDZZ = 0.0d0
!     do ilda = 1,inmx
!        do jlda = 1,inmx
!           write(200,*) ilda, jlda, DA_gpu(ilda,jlda), &
!                real( pft1(ilda,jlda)*conjg(pft2(ilda,jlda) ) ), &
!                DA_gpu(ilda,jlda) - real( pft1(ilda,jlda)*conjg(pft2(ilda,jlda) ) )
!           !sumDZZ = sumDZZ + DA_gpu(ilda,jlda) - real( pft1(ilda,jlda)*conjg(pft2(ilda,jlda) ) )
!        end do
!     end do

     r_gpu = sum(DA_gpu)
     !write(*,'(1x,a,1x,2f20.8)')"Full sum of the Double  precision matrix  DA_gpu is: ",r_gpu
!     write(*,'(1x,a,1x,f20.8)')"Full sum(DA_gpu-(Re(pft1*conjg(pft2)) of the Double precision matrix  sum is: ",sumDZZ

     r_gpu = calc_corr_dble(r_gpu,sumasq_gpu*sumbsq_gpu)
     write(*,*)
     write(*,*)"after loops r_gpu= ",r_gpu," sumasq_gpu= ",sumasq_gpu," sumbsq_gpu= ",sumbsq_gpu

     !write(*,'(i5,a,f10.5)')inmx,elps_corr

     !freeing the ressources on GPU
     err = cublas_free(devPtrA_DA)
     err = cublas_free(devPtrA_pft1)
     err = cublas_free(devPtrA_pft2)
     err = cublas_free(devPtrA_sumasq)
     err = cublas_free(devPtrA_sumbsq)

!Freeing the CPU memory 
     deallocate(pft1)  
     deallocate(pft2)
     deallocate(pft_gpu)
     deallocate(DA_gpu)


  end do

  !shutting down the environment
  call simple_cuda_shutdown()

  !shutting down the timers
  call stop_Alltimers_cpu()

  !end of greeting message
  call bye_gpu_magma()

end program testing_corr_gpu
