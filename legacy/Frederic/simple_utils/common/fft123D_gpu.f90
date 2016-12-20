!> \brief SIMPLE deviceQuery class for GPU
module fft123D_gpu
use, intrinsic :: iso_c_binding
use simple_defs
use simple_cuda_defs

implicit none
#define devptr_t integer*8

public :: cuFFT_gpu, hello_cuFFT_gpu, bye_cuFFT_gpu

!The c binding derived type
type, bind(c) :: cuFFT_devType
   integer(c_int) :: dim
   integer(c_int) :: npoints
end type cuFFT_devType

private 

type :: cuFFT_gpu

   private 
   type(cuFFT_devType) :: t_devFFT
   logical             :: existence_c_devFFT_gpu=.false. !< objects exist or not
   integer  :: nx
   ! double precision
   real(dp),allocatable :: x(:),y(:),z(:)
   real(dp),allocatable :: fx(:),fxy(:,:),fxyz(:,:,:)
   complex(dp),allocatable :: u(:),v(:),w(:)
   complex(dp),allocatable :: fu(:),fuv(:,:),fuvw(:,:,:)
   ! single precision
   real(sp),allocatable :: s_x(:),s_y(:),s_z(:)
   real(sp),allocatable :: s_fx(:),s_fxy(:,:),s_fxyz(:,:,:)
   complex(sp),allocatable :: s_u(:),s_v(:),s_w(:)
   complex(sp),allocatable :: s_fu(:),s_fuv(:,:),s_fuvw(:,:,:)

 contains
   !constructor
   procedure :: new_cuFFT_1D_ZZ_gpu   !Double precision
   procedure :: new_cuFFT_1D_DZ_gpu
   procedure :: new_cuFFT_2D_ZZ_gpu
   procedure :: new_cuFFT_2D_DZ_gpu
   procedure :: new_cuFFT_3D_ZZ_gpu
   procedure :: new_cuFFT_3D_DZ_gpu

   procedure :: new_cuFFT_3D_CC_gpu   !single precision
   procedure :: new_cuFFT_3D_SC_gpu
   procedure :: new_cuFFT_2D_CC_gpu
   procedure :: new_cuFFT_2D_SC_gpu

   !destructor
   procedure :: kill_cuFFT_gpu
   !Terminators
   procedure :: terminate
   !worker methods
   procedure :: gather_cuFFT
   procedure :: gather_fft_1D_Z2Z_gpu !double precision
   procedure :: gather_fft_2D_Z2Z_gpu
   procedure :: gather_fft_2D_D2Z_gpu
   procedure :: gather_fft_2D_Z2D_gpu
   procedure :: gather_fft_3D_Z2Z_gpu
   procedure :: gather_fft_3D_D2Z_gpu
   procedure :: gather_fft_3D_Z2D_gpu

   procedure :: gather_fft_2D_C2C_gpu !single precision
   procedure :: gather_fft_2D_S2C_gpu
   procedure :: gather_fft_2D_C2S_gpu
   procedure :: gather_fft_3D_C2C_gpu
   procedure :: gather_fft_3D_S2C_gpu
   procedure :: gather_fft_3D_C2S_gpu

   !setters
   procedure :: set_Ddata_1D
   procedure :: set_Zdata_1D
   procedure :: set_Ddata_2D
   procedure :: set_Zdata_2D
   procedure :: set_Ddata_3D
   procedure :: set_Zdata_3D
   !getters
   procedure :: get_Ddata_1D
   procedure :: get_Zdata_1D
   procedure :: get_x_1D
   procedure :: get_u_1D
   !initialisers
   procedure :: initialise_cuFFT
end type cuFFT_gpu

interface cuFFT_gpu
   module procedure constructor_cuFFT_gpu
end interface cuFFT_gpu

interface
#if defined (LINUX)
   !GPU
   ! 1D
   function gather_fft_1D_Z2Z_gpu_cpp(nx,fu,fh,sign)
     use simple_defs
     integer  :: nx
     integer  :: sign
     complex(dp) :: fu(*)
     complex(dp) :: fh(*)
     integer :: gather_fft_1D_Z2Z_gpu_cpp
   end function gather_fft_1D_Z2Z_gpu_cpp
   ! 2D
   function gather_fft_2D_Z2Z_gpu_cpp(nx,ny,fuv,fhp,sign)
     use simple_defs
     integer  :: nx,ny
     integer  :: sign
     complex(dp) :: fuv(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_Z2Z_gpu_cpp
   end function gather_fft_2D_Z2Z_gpu_cpp
   function gather_fft_2D_C2C_gpu_cpp(nx,ny,s_fuv,s_fhp,sign)
     use simple_defs
     integer  :: nx,ny
     integer  :: sign
     complex(sp) :: s_fuv(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_C2C_gpu_cpp
   end function gather_fft_2D_C2C_gpu_cpp
   function gather_fft_2D_D2Z_gpu_cpp(nx,ny,fxy,fhp)
     use simple_defs
     integer  :: nx,ny
     real(dp) :: fxy(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_D2Z_gpu_cpp
   end function gather_fft_2D_D2Z_gpu_cpp
   function gather_fft_2D_S2C_gpu_cpp(nx,ny,s_fxy,s_fhp)
     use simple_defs
     integer  :: nx,ny
     real(sp) :: s_fxy(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_S2C_gpu_cpp
   end function gather_fft_2D_S2C_gpu_cpp
   function gather_fft_2D_Z2D_gpu_cpp(nx,ny,fhp,fxy)
     use simple_defs
     integer  :: nx,ny
     real(dp) :: fxy(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_Z2D_gpu_cpp
   end function gather_fft_2D_Z2D_gpu_cpp
   function gather_fft_2D_C2S_gpu_cpp(nx,ny,s_fhp,s_fxy)
     use simple_defs
     integer  :: nx,ny
     real(sp) :: s_fxy(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_C2S_gpu_cpp
   end function gather_fft_2D_C2S_gpu_cpp
   ! 3D
   function gather_fft_3D_Z2Z_gpu_cpp(nx,ny,nz,fuvw,fhpq,sign)
     use simple_defs
     integer  :: nx,ny,nz
     integer  :: sign
     complex(dp) :: fuvw(nx,ny,*)
     complex(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_Z2Z_gpu_cpp
   end function gather_fft_3D_Z2Z_gpu_cpp
   function gather_fft_3D_C2C_gpu_cpp(nx,ny,nz,s_fuvw,s_fhpq,sign)
     use simple_defs
     integer  :: nx,ny,nz
     integer  :: sign
     complex(sp) :: s_fuvw(nx,ny,*)
     complex(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_C2C_gpu_cpp
   end function gather_fft_3D_C2C_gpu_cpp
   function gather_fft_3D_D2Z_gpu_cpp(nx,ny,nz,fxyz,fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     real(dp) :: fxyz(nx,ny,*)
     complex(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_D2Z_gpu_cpp
   end function gather_fft_3D_D2Z_gpu_cpp
   function gather_fft_3D_S2C_gpu_cpp(nx,ny,nz,s_fxyz,s_fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     real(sp) :: s_fxyz(nx,ny,*)
     complex(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_S2C_gpu_cpp
   end function gather_fft_3D_S2C_gpu_cpp
   function gather_fft_3D_Z2D_gpu_cpp(nx,ny,nz,fuvw,fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     complex(dp) :: fuvw(nx,ny,*)
     real(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_Z2D_gpu_cpp
   end function gather_fft_3D_Z2D_gpu_cpp
   function gather_fft_3D_C2S_gpu_cpp(nx,ny,nz,s_fuvw,s_fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     complex(sp) :: s_fuvw(nx,ny,*)
     real(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_C2S_gpu_cpp
   end function gather_fft_3D_C2S_gpu_cpp

#endif
end interface

contains 
  !CONSTRUCTORS
  !> \brief is a fft123D constructor
  function constructor_cuFFT_gpu(nx,ny,nz,fft_in) result(c_devFFT_gpu)
    type(cuFFT_gpu) :: c_devFFT_gpu
    integer,intent(inout),optional :: nx,ny,nz
    integer,intent(in) :: fft_in

    if ( present(nx) ) then
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_gpu%new_cuFFT_1D_DZ_gpu(nx)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_gpu%new_cuFFT_1D_ZZ_gpu(nx)
       end if
    else if ( present(nx) .and. present(ny) ) then
       !single precision
       if ( fft_in == FFT_S2C .or. fft_in == FFT_C2S ) then
          call c_devFFT_gpu%new_cuFFT_2D_SC_gpu(nx,ny)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_gpu%new_cuFFT_2D_CC_gpu(nx,ny)
       end if
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_gpu%new_cuFFT_2D_DZ_gpu(nx,ny)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_gpu%new_cuFFT_2D_ZZ_gpu(nx,ny)
       end if
    else if ( present(nx) .and. present(ny) .and. present(nz) ) then
       !single precision
       if ( fft_in == FFT_S2C .or. fft_in == FFT_C2S ) then
          call c_devFFT_gpu%new_cuFFT_3D_SC_gpu(nx,ny,nz)
       else if ( fft_in == FFT_C2C ) then
          call c_devFFT_gpu%new_cuFFT_3D_CC_gpu(nx,ny,nz)
       end if
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_gpu%new_cuFFT_3D_DZ_gpu(nx,ny,nz)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_gpu%new_cuFFT_3D_ZZ_gpu(nx,ny,nz)
       end if
    end if
   
  end function constructor_cuFFT_gpu
  ! 3D constructor method

  subroutine new_cuFFT_3D_CC_gpu(c_devFFT_gpu,nx,ny,nz)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny,nz

    !complex(sp)
    allocate( c_devFFT_gpu%s_u(nx))
    allocate( c_devFFT_gpu%s_v(ny))
    allocate( c_devFFT_gpu%s_w(nz))
    allocate( c_devFFT_gpu%s_fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_3D_CC_gpu
  subroutine new_cuFFT_3D_SC_gpu(c_devFFT_gpu,nx,ny,nz)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny,nz

    !real(sp)
    allocate( c_devFFT_gpu%s_x(nx))
    allocate( c_devFFT_gpu%s_y(ny))
    allocate( c_devFFT_gpu%s_z(nz))
    allocate( c_devFFT_gpu%fxyz(nx,ny,nz))
    !complex(sp)
    allocate( c_devFFT_gpu%s_u(nx))
    allocate( c_devFFT_gpu%s_v(ny))
    allocate( c_devFFT_gpu%s_w(nz))
    allocate( c_devFFT_gpu%s_fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_3D_SC_gpu
  subroutine new_cuFFT_3D_ZZ_gpu(c_devFFT_gpu,nx,ny,nz)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny,nz

    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%v(ny))
    allocate( c_devFFT_gpu%w(nz))
    allocate( c_devFFT_gpu%fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_3D_ZZ_gpu
  subroutine new_cuFFT_3D_DZ_gpu(c_devFFT_gpu,nx,ny,nz)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny,nz

    !real(dp)
    allocate( c_devFFT_gpu%x(nx))
    allocate( c_devFFT_gpu%y(ny))
    allocate( c_devFFT_gpu%z(nz))
    allocate( c_devFFT_gpu%fxyz(nx,ny,nz))
    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%v(ny))
    allocate( c_devFFT_gpu%w(nz))
    allocate( c_devFFT_gpu%fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_3D_DZ_gpu

  ! 2D constructor method 

  subroutine new_cuFFT_2D_CC_gpu(c_devFFT_gpu,nx,ny)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_gpu%s_u(nx))
    allocate( c_devFFT_gpu%s_v(ny))
    allocate( c_devFFT_gpu%s_fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_2D_CC_gpu
  subroutine new_cuFFT_2D_SC_gpu(c_devFFT_gpu,nx,ny)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices

    !real(dp)
    allocate( c_devFFT_gpu%s_x(nx))
    allocate( c_devFFT_gpu%s_y(ny))
    allocate( c_devFFT_gpu%s_fxy(nx,ny))
    !complex(dp)
    allocate( c_devFFT_gpu%s_u(nx))
    allocate( c_devFFT_gpu%s_v(ny))
    allocate( c_devFFT_gpu%s_fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_2D_SC_gpu
  subroutine new_cuFFT_2D_ZZ_gpu(c_devFFT_gpu,nx,ny)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%v(ny))
    allocate( c_devFFT_gpu%fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_2D_ZZ_gpu
  subroutine new_cuFFT_2D_DZ_gpu(c_devFFT_gpu,nx,ny)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices
    
    !real(dp)
    allocate( c_devFFT_gpu%x(nx))
    allocate( c_devFFT_gpu%y(ny))
    allocate( c_devFFT_gpu%fxy(nx,ny))
    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%v(ny))
    allocate( c_devFFT_gpu%fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_2D_DZ_gpu

  ! 1D constructor method 
  subroutine new_cuFFT_1D_ZZ_gpu(c_devFFT_gpu,nx)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%fu(nx))

    call gather_cuFFT(c_devFFT_gpu,t_devFFT)

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_1D_ZZ_gpu
  subroutine new_cuFFT_1D_DZ_gpu(c_devFFT_gpu,nx)
    class(cuFFT_gpu) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: nx
    !kill the existing object before allocating a new one
    call c_devFFT_gpu%kill_cuFFT_gpu

    !allocating arrays for details for n devices
    
    !real(dp)
    allocate( c_devFFT_gpu%x(nx))
    allocate( c_devFFT_gpu%fx(nx))
    !complex(dp)
    allocate( c_devFFT_gpu%u(nx))
    allocate( c_devFFT_gpu%fu(nx))

    call gather_cuFFT(c_devFFT_gpu,t_devFFT)

    !set the existence of the object to true
    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine new_cuFFT_1D_DZ_gpu

  !WORKERS METHODS
  subroutine gather_cuFFT(c_devFFT_gpu,t_devFFT)
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    !TODO: general gather method to parse general info may be removed later
    return
  end subroutine gather_cuFFT

  ! 3D data
  ! getting the fourier transform on GPU for
  ! Z2Z
  subroutine gather_fft_3D_Z2Z_gpu(c_devFFT_gpu,nx,ny,nz,fuvw,fhpq,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz,sign
    complex(dp) :: fuvw(nx,ny,*)
    complex(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_Z2Z_gpu_c
       function gather_fft_3D_Z2Z_gpu_c(nx,ny,nz,fuvw,fhpq,sign)
         use simple_defs
         integer  :: nx,ny,nz
         integer  :: sign
         complex(dp) :: fuvw(nx,ny,*)
         complex(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_Z2Z_gpu_c
       end function gather_fft_3D_Z2Z_gpu_c
    end interface external_c_function_gather_fft_3D_Z2Z_gpu_c
    
#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_Z2Z_gpu_c(nx,ny,nz,fuvw,fhpq,sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_Z2Z_gpu_cpp(nx,ny,nz,fuvw,fhpq,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_3D_Z2Z_gpu
  ! C2C
  subroutine gather_fft_3D_C2C_gpu(c_devFFT_gpu,nx,ny,nz,s_fuvw,s_fhpq,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz,sign
    complex(sp) :: s_fuvw(nx,ny,*)
    complex(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_C2C_gpu_c
       function gather_fft_3D_C2C_gpu_c(nx,ny,nz,s_fuvw,s_fhpq,sign)
         use simple_defs
         integer  :: nx,ny,nz
         integer  :: sign
         complex(sp) :: s_fuvw(nx,ny,*)
         complex(sp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_C2C_gpu_c
       end function gather_fft_3D_C2C_gpu_c
    end interface external_c_function_gather_fft_3D_C2C_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_C2C_gpu_c(nx,ny,nz,s_fuvw,s_fhpq,sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_C2C_gpu_cpp(nx,ny,nz,s_fuvw,s_fhpq,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_3D_C2C_gpu
  ! D2Z
  subroutine gather_fft_3D_D2Z_gpu(c_devFFT_gpu,nx,ny,nz,fxyz,fhpq)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz
    real(dp) :: fxyz(nx,ny,*)
    complex(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_D2Z_gpu_c
       function gather_fft_3D_D2Z_gpu_c(nx,ny,nz,fxyz,fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         real(dp) :: fxyz(nx,ny,*)
         complex(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_D2Z_gpu_c
       end function gather_fft_3D_D2Z_gpu_c
    end interface external_c_function_gather_fft_3D_D2Z_gpu_c
#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_D2Z_gpu_c(nx,ny,nz,fxyz,fhpq)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_D2Z_gpu_cpp(nx,ny,nz,fxyz,fhpq)
#else
    call c_devFFT_gpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_D2Z_gpu
  ! S2C
  subroutine gather_fft_3D_S2C_gpu(c_devFFT_gpu,nx,ny,nz,s_fxyz,s_fhpq)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz
    real(sp) :: s_fxyz(nx,ny,*)
    complex(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_S2C_gpu_c
       function gather_fft_3D_S2C_gpu_c(nx,ny,nz,s_fxyz,s_fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         real(dp) :: s_fxyz(nx,ny,*)
         complex(dp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_S2C_gpu_c
       end function gather_fft_3D_S2C_gpu_c
    end interface external_c_function_gather_fft_3D_S2C_gpu_c
#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_S2C_gpu_c(nx,ny,nz,s_fxyz,s_fhpq)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_S2C_gpu_cpp(nx,ny,nz,s_fxyz,s_fhpq)
#else
    call c_devFFT_gpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_S2C_gpu
  ! Z2D
  subroutine gather_fft_3D_Z2D_gpu(c_devFFT_gpu,nx,ny,nz,fuvw,fhpq)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz
    complex(dp) :: fuvw(nx,ny,*)
    real(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_Z2D_gpu_c
       function gather_fft_3D_Z2D_gpu_c(nx,ny,nz,fuvw,fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         complex(dp) :: fuvw(nx,ny,*)
         real(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_Z2D_gpu_c
       end function gather_fft_3D_Z2D_gpu_c
    end interface external_c_function_gather_fft_3D_Z2D_gpu_c
#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_Z2D_gpu_c(nx,ny,nz,fuvw,fhpq)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_Z2D_gpu_cpp(nx,ny,nz,fuvw,fhpq)
#else
    call c_devFFT_gpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_Z2D_gpu
  ! C2S
  subroutine gather_fft_3D_C2S_gpu(c_devFFT_gpu,nx,ny,nz,s_fuvw,s_fhpq)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,nz
    complex(sp) :: s_fuvw(nx,ny,*)
    real(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_C2S_gpu_c
       function gather_fft_3D_C2S_gpu_c(nx,ny,nz,s_fuvw,s_fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         complex(sp) :: s_fuvw(nx,ny,*)
         real(sp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_C2S_gpu_c
       end function gather_fft_3D_C2S_gpu_c
    end interface external_c_function_gather_fft_3D_C2S_gpu_c
#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_3D_C2S_gpu_c(nx,ny,nz,s_fuvw,s_fhpq)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_3D_C2S_gpu_cpp(nx,ny,nz,s_fuvw,s_fhpq)
#else
    call c_devFFT_gpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_C2S_gpu

  ! 2D data
  !getting the fourier transform on GPU for
  ! Z2Z
  subroutine gather_fft_2D_Z2Z_gpu(c_devFFT_gpu,nx,ny,fuv,fhp,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,sign
    complex(dp) :: fuv(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_Z2Z_gpu_c
       function gather_fft_2D_Z2Z_gpu_c(nx, ny,fuv,fhp,sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         complex(dp) :: fuv(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_Z2Z_gpu_c
       end function gather_fft_2D_Z2Z_gpu_c
    end interface external_c_function_gather_fft_2D_Z2Z_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_Z2Z_gpu_c(nx, ny, fuv, fhp, sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_Z2Z_gpu_cpp(nx, ny,fuv,fhp,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_Z2Z_gpu
  ! C2C
  subroutine gather_fft_2D_C2C_gpu(c_devFFT_gpu,nx,ny,s_fuv,s_fhp,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny,sign
    complex(sp) :: s_fuv(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_C2C_gpu_c
       function gather_fft_2D_C2C_gpu_c(nx, ny,s_fuv,s_fhp,sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         complex(sp) :: s_fuv(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_C2C_gpu_c
       end function gather_fft_2D_C2C_gpu_c
    end interface external_c_function_gather_fft_2D_C2C_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_C2C_gpu_c(nx, ny, s_fuv, s_fhp, sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_C2C_gpu_cpp(nx, ny,s_fuv,s_fhp,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_C2C_gpu
  ! D2Z
  subroutine gather_fft_2D_D2Z_gpu(c_devFFT_gpu,nx,ny,fxy,fhp)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny
    real(dp) :: fxy(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_D2Z_gpu_c
       function gather_fft_2D_D2Z_gpu_c(nx, ny,fxy,fhp)
         use simple_defs
         integer  :: nx,ny
         real(dp) :: fxy(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_D2Z_gpu_c
       end function gather_fft_2D_D2Z_gpu_c
    end interface external_c_function_gather_fft_2D_D2Z_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_D2Z_gpu_c(nx, ny, fxy, fhp)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_D2Z_gpu_cpp(nx, ny,fxy,fhp)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_D2Z_gpu
  ! S2C
  subroutine gather_fft_2D_S2C_gpu(c_devFFT_gpu,nx,ny,s_fxy,s_fhp)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny
    real(sp) :: s_fxy(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_S2C_gpu_c
       function gather_fft_2D_S2C_gpu_c(nx, ny,s_fxy,s_fhp)
         use simple_defs
         integer  :: nx,ny
         real(sp) :: s_fxy(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_S2C_gpu_c
       end function gather_fft_2D_S2C_gpu_c
    end interface external_c_function_gather_fft_2D_S2C_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_S2C_gpu_c(nx, ny, s_fxy, s_fhp)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_S2C_gpu_cpp(nx, ny,s_fxy,s_fhp)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_S2C_gpu
  ! Z2D
  subroutine gather_fft_2D_Z2D_gpu(c_devFFT_gpu,nx,ny,fhp,fxy)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny
    real(dp) :: fxy(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_Z2D_gpu_c
       function gather_fft_2D_Z2D_gpu_c(nx, ny,fhp,fxy)
         use simple_defs
         integer  :: nx,ny
         real(dp) :: fxy(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_Z2D_gpu_c
       end function gather_fft_2D_Z2D_gpu_c
    end interface external_c_function_gather_fft_2D_Z2D_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_Z2D_gpu_c(nx, ny, fhp, fxy)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_Z2D_gpu_cpp(nx, ny,fhp,fxy)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_Z2D_gpu
  ! C2S
  subroutine gather_fft_2D_C2S_gpu(c_devFFT_gpu,nx,ny,s_fhp,s_fxy)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,ny
    real(sp) :: s_fxy(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_C2S_gpu_c
       function gather_fft_2D_C2S_gpu_c(nx, ny,s_fhp,s_fxy)
         use simple_defs
         integer  :: nx,ny
         real(sp) :: s_fxy(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_C2S_gpu_c
       end function gather_fft_2D_C2S_gpu_c
    end interface external_c_function_gather_fft_2D_C2S_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_2D_C2S_gpu_c(nx, ny, s_fhp, s_fxy)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_2D_C2S_gpu_cpp(nx, ny,s_fhp,s_fxy)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_C2S_gpu

  ! 1D data
  !getting the fourier transform on GPU for
  subroutine gather_fft_1D_Z2Z_gpu(c_devFFT_gpu,nx,fu,fh,sign)
    use simple_defs
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer  :: nx,sign
    complex(dp) :: fu(*)
    complex(dp) :: fh(*)
    integer :: rc !return code
    interface external_c_function_gather_fft_1D_Z2Z_gpu_c
       function gather_fft_1D_Z2Z_gpu_c(nx,fu,fh,sign)
         use simple_defs
         integer  :: nx
         integer  :: sign
         complex(dp) :: fu(*)
         complex(dp) :: fh(*)
         integer :: gather_fft_1D_Z2Z_gpu_c
       end function gather_fft_1D_Z2Z_gpu_c
    end interface external_c_function_gather_fft_1D_Z2Z_gpu_c

#if defined (MACOSX) && defined (CUDA)
    rc = gather_fft_1D_Z2Z_gpu_c(nx,fu,fh,sign)
#elif defined (LINUX) && defined (CUDA)
    rc = gather_fft_1D_Z2Z_gpu_cpp(nx,fu,fh,sign)
#else
    call c_devFFT_gpu%terminate()
#endif

    return
  end subroutine gather_fft_1D_Z2Z_gpu

  !SETTERS
  ! 1D
  subroutine set_Ddata_1D(c_devFFT_gpu,nx,x,Ddata1D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer  :: nx
    real(dp) :: x(*)
    real(dp) :: Ddata1D(*)
    c_devFFT_gpu%x(1:nx) = x(1:nx)
    c_devFFT_gpu%fx(1:nx) = Ddata1D(1:nx)
    return
  end subroutine set_Ddata_1D
  subroutine set_Zdata_1D(c_devFFT_gpu,nx,u,Zdata1D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer  :: nx
    complex(dp) :: u(*)
    complex(dp) :: Zdata1D(*)
    c_devFFT_gpu%u(1:nx) = u(1:nx)
    c_devFFT_gpu%fu(1:nx) = Zdata1D(1:nx)
    return
  end subroutine set_Zdata_1D
  ! 2D
  subroutine set_Ddata_2D(c_devFFT_gpu,nx,ny,Ddata2D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer :: nx,ny
    real(dp) :: Ddata2D(nx,*)
    c_devFFT_gpu%fxy(1:nx,1:ny) = Ddata2D(1:nx,1:ny)
    return
  end subroutine set_Ddata_2D
  subroutine set_Zdata_2D(c_devFFT_gpu,nx,ny,Zdata2D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer :: nx,ny
    complex(dp) :: Zdata2D(nx,*)
    c_devFFT_gpu%fuv(1:nx,1:ny) = Zdata2D(1:nx,1:ny)
    return
  end subroutine set_Zdata_2D
  ! 3D
  subroutine set_Ddata_3D(c_devFFT_gpu,nx,ny,nz,Ddata3D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer :: nx,ny,nz
    real(dp) :: Ddata3D(nx,ny,*)
    c_devFFT_gpu%fxyz(1:nx,1:ny,1:nz) = Ddata3D(1:nx,1:ny,1:nz)
    return
  end subroutine set_Ddata_3D
  subroutine set_Zdata_3D(c_devFFT_gpu,nx,ny,nz,Zdata3D)
    implicit none
    class(cuFFT_gpu),intent(inout) :: c_devFFT_gpu
    integer :: nx,ny,nz
    complex(dp) :: Zdata3D(nx,ny,*)
    c_devFFT_gpu%fuvw(1:nx,1:ny,1:nz) = Zdata3D(1:nx,1:ny,1:nz)
    return
  end subroutine set_Zdata_3D

  !GETTERS
  !double precision
  function get_Ddata_1D(c_devFFT_gpu,nx) result(Ddata1D_out)
    class(cuFFT_gpu),intent(in) :: c_devFFT_gpu
    integer :: nx
    real(dp) :: Ddata1D_out(1:nx)
    Ddata1D_out(1:nx) = c_devFFT_gpu%fx(1:nx)
  end function get_Ddata_1D
  function get_x_1D(c_devFFT_gpu,nx) result(x_out)
    class(cuFFT_gpu),intent(in) :: c_devFFT_gpu
    integer :: nx
    real(dp) :: x_out(1:nx)
    x_out(1:nx) = c_devFFT_gpu%x(1:nx)
  end function get_x_1D
  ! complex double precision
  function get_Zdata_1D(c_devFFT_gpu,nx) result(Zdata1D_out)
    class(cuFFT_gpu),intent(in) :: c_devFFT_gpu
    integer :: nx
    complex(dp) :: Zdata1D_out(1:nx)
    Zdata1D_out(1:nx) = c_devFFT_gpu%fu(1:nx)
  end function get_Zdata_1D
  function get_u_1D(c_devFFT_gpu,nx) result(u_out)
    class(cuFFT_gpu),intent(in) :: c_devFFT_gpu
    integer :: nx
    complex(dp) :: u_out(1:nx)
    u_out(1:nx) = c_devFFT_gpu%u(1:nx)
  end function get_u_1D

  !INITIALISERS
  subroutine initialise_cuFFT(c_devFFT_gpu,idev,t_devFFT)
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    type(cuFFT_devType) :: t_devFFT
    integer :: idev
    interface external_c_function_init
    end interface external_c_function_init

    idev = 0
    !TODO: add setters for nx,ny,nz dim

#if defined (MACOSX) && defined (CUDA)
#elif defined (LINUX) && defined (CUDA)
#else
    call c_devFFT_gpu%terminate()
#endif

    c_devFFT_gpu%existence_c_devFFT_gpu = .true.

    return
  end subroutine initialise_cuFFT

  !GREETERS
  !> \brief hello greeting routine for the object
  subroutine hello_cuFFT_gpu(err)
    implicit none
    integer :: err
    !start of the execution commands
    write(*,*) "Hello cuFFT GPU world"
    write(*,*)
    err = 0
    return
  end subroutine hello_cuFFT_gpu

  !> \brief bye greeting routine for the object
  subroutine bye_cuFFT_gpu()
    implicit none
    !start of the execution commands
    write(*,*) "Bye cuFFT GPU world"
    write(*,*)
    return
  end subroutine bye_cuFFT_gpu

  !TERMINATORS because of incompatiblilties errors
  subroutine terminate(c_devFFT_gpu)
    class(cuFFT_gpu) :: c_devFFT_gpu
    call c_devFFT_gpu%kill_cuFFT_gpu
    write(*,*)"**************************WARNING*******************************"
    write(*,*)"*There are no GPU devices available for computation            *"
    write(*,*)"*You need to check that your hardware is suitable and has GPU  *"
    write(*,*)"*computational capacities or that you have compiled with       *"
    write(*,*)"*-DCUDA to acces the CUDA environment computation.             *"
    write(*,*)"****************************************************************"
    call bye_cuFFT_gpu()
    stop
    return
  end subroutine terminate

  !DESTRUCTOR
  !> \brief is a cuFFT_gpu destructor
  subroutine kill_cuFFT_gpu(c_devFFT_gpu)
    class(cuFFT_gpu), intent(inout) :: c_devFFT_gpu
    if ( c_devFFT_gpu%existence_c_devFFT_gpu) then

       if(allocated(c_devFFT_gpu%x )) deallocate(c_devFFT_gpu%x )
       if(allocated(c_devFFT_gpu%y )) deallocate(c_devFFT_gpu%y )
       if(allocated(c_devFFT_gpu%z )) deallocate(c_devFFT_gpu%z )

       if(allocated(c_devFFT_gpu%s_x )) deallocate(c_devFFT_gpu%s_x )
       if(allocated(c_devFFT_gpu%s_y )) deallocate(c_devFFT_gpu%s_y )
       if(allocated(c_devFFT_gpu%s_z )) deallocate(c_devFFT_gpu%s_z )

       if(allocated(c_devFFT_gpu%fx)) deallocate(c_devFFT_gpu%fx)
       if(allocated(c_devFFT_gpu%fxy)) deallocate(c_devFFT_gpu%fxy)
       if(allocated(c_devFFT_gpu%fxyz)) deallocate(c_devFFT_gpu%fxyz)

       if(allocated(c_devFFT_gpu%s_fx)) deallocate(c_devFFT_gpu%s_fx)
       if(allocated(c_devFFT_gpu%s_fxy)) deallocate(c_devFFT_gpu%s_fxy)
       if(allocated(c_devFFT_gpu%s_fxyz)) deallocate(c_devFFT_gpu%s_fxyz)

       if(allocated(c_devFFT_gpu%u )) deallocate(c_devFFT_gpu%u )
       if(allocated(c_devFFT_gpu%v )) deallocate(c_devFFT_gpu%v )
       if(allocated(c_devFFT_gpu%w )) deallocate(c_devFFT_gpu%w )

       if(allocated(c_devFFT_gpu%s_u )) deallocate(c_devFFT_gpu%s_u )
       if(allocated(c_devFFT_gpu%s_v )) deallocate(c_devFFT_gpu%s_v )
       if(allocated(c_devFFT_gpu%s_w )) deallocate(c_devFFT_gpu%s_w )

       if(allocated(c_devFFT_gpu%fu)) deallocate(c_devFFT_gpu%fu)
       if(allocated(c_devFFT_gpu%fuv)) deallocate(c_devFFT_gpu%fuv)
       if(allocated(c_devFFT_gpu%fuvw)) deallocate(c_devFFT_gpu%fuvw)

       if(allocated(c_devFFT_gpu%s_fu)) deallocate(c_devFFT_gpu%s_fu)
       if(allocated(c_devFFT_gpu%s_fuv)) deallocate(c_devFFT_gpu%s_fuv)
       if(allocated(c_devFFT_gpu%s_fuvw)) deallocate(c_devFFT_gpu%s_fuvw)

       c_devFFT_gpu%existence_c_devFFT_gpu = .false.

    end if
    return
  end subroutine kill_cuFFT_gpu

end module fft123D_gpu
