!>  \brief SIMPLE polarft class for GPU
module simple_polarft_gpu
use simple_defs
use simple_cuda_defs !TODO: fixe the cuda linkage
use simple_polarft
use simple_jiffys, only: alloc_err
!use simple_jiffys_gpu, only: alloc_gpu_err :TODO: to be implemented
implicit none
#define devptr_t integer*8

public :: polarft_gpu
private
    integer               :: ldim(3)=0 !< logical dimensions of original cartesian image
type, extends(polarft) :: polarft_gpu
   private
   !TODO: add the private variables
   logical               :: existence_gpu=.false. !< objects exist or not
   integer               :: nradial_gpu   !< # radial vectors (angle index)(=# components in each shell))
   integer               :: khp_gpu       !< high-pass frequency limit 
   integer               :: klp_gpu       !< low-pass frequency limit
   integer               :: ring2_gpu     !< radius of molecule
   complex, allocatable  :: pft_gpu(:,:)  !< first coord is angle (1:nradia

 contains
   !TODO: add the methods to implement on GPU
   !constructor
   procedure :: new_gpu
   !destructor
   procedure :: kill_gpu
end type polarft_gpu

interface polarft_gpu
   module procedure constructor_gpu
end interface polarft_gpu
! CLASS PARAMETERS/VARIABLES
real, allocatable  :: polar_angtab(:)     !< angles (in degrees)
real, allocatable  :: polar_coords(:,:,:) !< polar coordinates (Cartesian)
contains

  ! CONSTRUCTORS

  !> \brief is a polar Fourier constructor
  function constructor_gpu( kfromto, ring2 ) result(self_gpu)
    integer, intent(in) :: kfromto(2), ring2
    type(polarft_gpu)  :: self_gpu
    !implement the constructor
    call self_gpu%new_gpu(kfromto, ring2)
  end function constructor_gpu
  
  !>  \brief  is a polar Fourier transform constructor
  subroutine new_gpu( self_gpu, kfromto, ring2 )
    use simple_math_gpu, only: gen_polar_coords_gpu
    class(polarft_gpu), intent(inout) :: self_gpu
    integer, intent(in)           :: kfromto(2), ring2
    integer                       :: alloc_stat
    call self_gpu%kill_gpu
    if( kfromto(2)-kfromto(1) < 1 )then
       write(*,*) 'fromto:', kfromto(1), kfromto(2)
       stop 'resolution range to narrow; new; simple_polarft'
    endif
    if( ring2 > 0 )then
    else
       stop 'ring2 must be > 0; new; simple_polarft'
    endif
    call gen_polar_coords_gpu(kfromto, ring2, polar_coords, polar_angtab)

    self_gpu%khp_gpu     = kfromto(1)
    self_gpu%klp_gpu     = kfromto(2)
    self_gpu%ring2_gpu   = ring2
    self_gpu%nradial_gpu = size(polar_angtab)

    if( self_gpu%nradial_gpu > 0 )then
    else
       stop 'ring2 must be > 0; new; simple_polarft'
    endif
    allocate(self_gpu%pft_gpu(self_gpu%nradial_gpu,self_gpu%khp_gpu:self_gpu%klp_gpu), stat=alloc_stat)
    call alloc_err("In: new; polarft", alloc_stat)
    self_gpu%existence_gpu = .true.
    return
  end subroutine new_gpu
    
  !>  \brief  is for correlating two polar images
  !!          parameterization is over rotations
  function corr_gpu( self1, rot1, self2, rot2 ) result( r )
    use simple_math, only: calc_corr_dble
    class(polarft), intent(in) :: self1
    class(polarft), intent(in) :: self2
    integer, intent(in)        :: rot1, rot2
    real                       :: r, sumasq, sumbsq
    integer                    :: i1, i2
    !local variables
    integer                    :: n_pft
    integer                    :: m,n,k
    integer                    :: lda, ldb, ldc
    real(dp)                   :: dble_alpha,dble_beta
    !gpu variables for the calculation of the r, sumasq and sumbsq

    real(dp)                   :: r_gpu, sumasq_gpu, sumbsq_gpu
    real(dp), allocatable      :: DA_gpu(:,:)
    complex(dp), allocatable   :: pft1(:,:)          !< first coord is angle (1:nradia
    complex(dp), allocatable   :: pft2(:,:)          !< first coord is angle (1:nradia

    !the device pointers 
    devptr_t                   ::  devPtrA_D1
    devptr_t                   ::  devPtrA_D2
    devptr_t                   ::  devPtrA_D3
    devptr_t                   ::  devPtrA_pft1
    devptr_t                   ::  devPtrA_pft2

    !function calls
    integer, external          :: get_nradial
    integer, external          :: get_khp
    integer, external          :: get_klp
    integer, external          :: get_ring2

    !return code
    integer                    :: err
#if defined (CUDA)

    err = 0

    if( self1.eqdims.self2 )then
       i1     = rot1
       i2     = rot2
       sumasq = 0.
       sumbsq = 0.
       r      = 0.

       !***************GPU corr****************************************'

       !proceeding to GPU calculation
       dble_alpha = 1.0d0
       dble_beta = 0.0d0
       !allocating memory
       allocate(pft1(get_nradial(self1),get_khp(self1):get_klp(self1)))
       allocate(pft2(get_nradial(self2),get_khp(self2):get_klp(self2)))

       !n_pft = (self1%nradial/2) * (self1%klp - self1%khp)
       n_pft = ( get_nradial( self1 ) / 2 ) * ( get_klp(self1) - get_khp(self1) )
       m = get_nradial( self1 ) / 2         !TODO: (self1%nradial/2)       the size of the matrix
       n = get_klp(self1) - get_khp(self1)  !TODO: (self1%klp - self1%khp) the size of the matrices
       k = get_khp(self1)                   !TODO: (self1%klp)

       lda = m   !TODO: must fix to the correct values
       ldb = lda !TODO: must fix the lda 
       ldc = lda !TODO: must fix the lda

       err = cublas_alloc(m*n, size_of_double_complex, devPtrA_pft1)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
       err = cublas_alloc(m*n, size_of_double_complex, devPtrA_pft2)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       !setting up the pft 1 and 2 matrix on device
!       pft1 = dble(get_pft(self1))
!       err = cublas_set_matrix (m, n, size_of_double_complex, pft1, lda, devPtrA_pft1, lda )
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
!       pft2 = dble(get_pft(self2))
!       err = cublas_set_matrix (m, n, size_of_double_complex, pft2, ldb, devPtrA_pft2, ldb )
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       !Now calluating the r value
       allocate(DA_gpu(m,n))
       err = cublas_alloc(m*n, size_of_double, devPtrA_D1)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
       err = cublas_set_matrix (m, n, size_of_double, DA_gpu, lda, devPtrA_D1, lda )
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       call zz2dgemm_ElmtWs_tesla_gpu("N", "N", &
                                      m, n, k,  &
                                      dble_alpha,  &
                                      devPtrA_pft1, lda,  &
                                      devPtrA_pft2, ldb,  &
                                      dble_beta,  &
                                      devPtrA_D1, ldc)

       DA_gpu = 0.0
       err = cublas_get_matrix ( m, n, size_of_double, devPtrA_D1, lda, DA_gpu, lda)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
       
       r_gpu = sum(DA_gpu)
       !freein the ressources on device for the first r calculation.
       err = cublas_free(devPtrA_D1)
       deallocate(DA_gpu)

       !***********Now calculate the sumsq(a and b)*********************
       !******sumasq*****
       allocate(DA_gpu(m,n))

       err = cublas_alloc(m*n, size_of_double, devPtrA_D2)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)
       err = cublas_set_matrix (m, n, size_of_double, DA_gpu, lda, devPtrA_D2, lda )
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       call zz2dgemm_ElmtWs_tesla_sumsq_gpu("N", "N", &
                                            m, n, k,  &
                                            dble_alpha,  &
                                            devPtrA_pft1, lda,  &
                                            devPtrA_pft2, ldb,  &
                                            dble_beta,  &
                                            devPtrA_D2, ldc)

       DA_gpu = 0.0
       err = cublas_get_matrix ( m, n, size_of_double, devPtrA_D2, lda, DA_gpu, lda)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       sumasq_gpu = 0.0d0
       sumasq_gpu = sum(DA_gpu)
       err = cublas_free(devPtrA_D2)
       deallocate(DA_gpu)

       !******sumbsq*****
       allocate(DA_gpu(m,n))
       err = cublas_alloc(m*n, size_of_double, devPtrA_D3)
       err = cublas_set_matrix (m, n, size_of_double, DA_gpu, lda, devPtrA_D3, lda )

       call zz2dgemm_ElmtWs_tesla_sumsq_gpu("N", "N", &
                                            m, n, k,  &
                                            dble_alpha,  &
                                            devPtrA_pft1, lda,  &
                                            devPtrA_pft2, ldb,  &
                                            dble_beta,  &
                                            devPtrA_D3, ldc)

       DA_gpu = 0.0
       err = cublas_get_matrix ( m, n, size_of_double, devPtrA_D3, lda, DA_gpu, lda)
       if ( err .ne. 0 ) call simple_cudblas_stat_return(err)

       sumbsq_gpu = 0.0d0
       sumbsq_gpu = sum(DA_gpu)
       err = cublas_free(devPtrA_D3)

       deallocate(DA_gpu)

       !******Now calculating the correlator*****
       r = 0.0
       r = real(calc_corr_dble(r_gpu,sumasq_gpu*sumbsq_gpu))

       !***************************************************************

       !freeing the ressources on GPU
       err = cublas_free(devPtrA_pft1)
       err = cublas_free(devPtrA_pft2)

       !r = calc_corr(r,sumasq*sumbsq)
    else
       write(*,*) get_ring2(self1), get_khp(self1), get_klp(self1)
       write(*,*) get_ring2(self2), get_khp(self2), get_klp(self2)
       stop 'not equal dims; corr'
    endif

#else
    write(*,*)"**************************WARNING******************************"
    write(*,*)"You need to compile with -DCUDA                                "
    write(*,*)"to acces the CUDA environment computation using GPU            "
    write(*,*)"switching back to the CPU version of corr function             "
    write(*,*)"***************************************************************"

    r = self1%corr_slow( rot1, self2, rot2 )

#endif

    return
  end function corr_gpu

  !> \brief is a polar_gpu desturctor
  subroutine kill_gpu(self_gpu)
    class(polarft_gpu), intent(inout) :: self_gpu
    !TODO: implement the destructor

    if (self_gpu%existence_gpu ) then
       deallocate(self_gpu%pft_gpu)
       if( allocated(polar_angtab) ) deallocate(polar_angtab)
       if( allocated(polar_coords) ) deallocate(polar_coords)
       self_gpu%existence_gpu = .false.
    end if

    return
  end subroutine kill_gpu

end module simple_polarft_gpu
