!> \brief SIMPLE deviceQuery class for GPU
module fft123D_cpu
use, intrinsic :: iso_c_binding
use simple_defs
use simple_cuda_defs
use simple_fftw3

implicit none

public :: fft_cpu, hello_fft_cpu, bye_fft_cpu

!The c binding derived type
type, bind(c) :: fft_devType
   integer(c_int) :: dim
   integer(c_int) :: npoints
end type fft_devType

private 

type :: fft_cpu

   private 
   type(fft_devType) :: t_devFFT
   logical             :: existence_c_devFFT_cpu=.false. !< objects exist or not
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
   procedure :: new_fft_1D_ZZ_cpu ! double precision
   procedure :: new_fft_1D_DZ_cpu
   procedure :: new_fft_2D_ZZ_cpu
   procedure :: new_fft_2D_DZ_cpu
   procedure :: new_fft_3D_ZZ_cpu
   procedure :: new_fft_3D_DZ_cpu

   procedure :: new_fft_3D_CC_cpu !single precision
   procedure :: new_fft_3D_SC_cpu
   procedure :: new_fft_2D_CC_cpu
   procedure :: new_fft_2D_SC_cpu

   !destructor
   procedure :: kill_fft_cpu
   !Terminators
   procedure :: terminate
   !worker methods
   procedure :: gather_fft
   procedure :: gather_fft_1D_Z2Z_cpu ! double precision
   procedure :: gather_fft_2D_Z2Z_cpu
   procedure :: gather_fft_2D_D2Z_cpu
   procedure :: gather_fft_2D_Z2D_cpu
   procedure :: gather_fft_3D_Z2Z_cpu
   procedure :: gather_fft_3D_D2Z_cpu
   procedure :: gather_fft_3D_Z2D_cpu

   procedure :: gather_fft_3D_C2C_cpu ! single precision
   procedure :: gather_fft_3D_S2C_cpu
   procedure :: gather_fft_3D_C2S_cpu
   procedure :: gather_fft_2D_C2C_cpu
   procedure :: gather_fft_2D_S2C_cpu
   procedure :: gather_fft_2D_C2S_cpu

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
   procedure :: initialise_fft
end type fft_cpu

interface fft_cpu
   module procedure constructor_fft_cpu
end interface fft_cpu

interface

#if defined (LINUX)
   !CPU
   ! 1D
   function gather_fft_1D_Z2Z_cpu_cpp(nx,fu,fh,sign)
     use simple_defs
     integer  :: nx
     integer  :: sign
     complex(dp) :: fu(*)
     complex(dp) :: fh(*)
     integer :: gather_fft_1D_Z2Z_cpu_cpp
   end function gather_fft_1D_Z2Z_cpu_cpp

   ! 2D
   function gather_fft_2D_Z2Z_cpu_cpp(nx,ny,fuv,fhp,sign)
     use simple_defs
     integer  :: nx,ny
     integer  :: sign
     complex(dp) :: fuv(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_Z2Z_cpu_cpp
   end function gather_fft_2D_Z2Z_cpu_cpp
   function gather_fft_2D_C2C_cpu_cpp(nx,ny,s_fuv,s_fhp,sign)
     use simple_defs
     integer  :: nx,ny
     integer  :: sign
     complex(sp) :: s_fuv(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_C2C_cpu_cpp
   end function gather_fft_2D_C2C_cpu_cpp
   function gather_fft_2D_D2Z_cpu_cpp(nx,ny,fxy,fhp)
     use simple_defs
     integer  :: nx,ny
     real(dp) :: fxy(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_D2Z_cpu_cpp
   end function gather_fft_2D_D2Z_cpu_cpp
   function gather_fft_2D_S2C_cpu_cpp(nx,ny,s_fxy,s_fhp)
     use simple_defs
     integer  :: nx,ny
     real(sp) :: s_fxy(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_S2C_cpu_cpp
   end function gather_fft_2D_S2C_cpu_cpp
   function gather_fft_2D_Z2D_cpu_cpp(nx,ny,fhp,fxy)
     use simple_defs
     integer  :: nx,ny
     real(dp) :: fxy(nx,*)
     complex(dp) :: fhp(nx,*)
     integer :: gather_fft_2D_Z2D_cpu_cpp
   end function gather_fft_2D_Z2D_cpu_cpp
   function gather_fft_2D_C2S_cpu_cpp(nx,ny,s_fhp,s_fxy)
     use simple_defs
     integer  :: nx,ny
     real(sp) :: s_fxy(nx,*)
     complex(sp) :: s_fhp(nx,*)
     integer :: gather_fft_2D_C2S_cpu_cpp
   end function gather_fft_2D_C2S_cpu_cpp
   ! 3D
   function gather_fft_3D_Z2Z_cpu_cpp(nx,ny,nz,fuvw,fhpq,sign)
     use simple_defs
     integer  :: nx,ny,nz
     integer  :: sign
     complex(dp) :: fuvw(nx,ny,*)
     complex(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_Z2Z_cpu_cpp
   end function gather_fft_3D_Z2Z_cpu_cpp
   function gather_fft_3D_C2C_cpu_cpp(nx,ny,nz,s_fuvw,s_fhpq,sign)
     use simple_defs
     integer  :: nx,ny,nz
     integer  :: sign
     complex(sp) :: s_fuvw(nx,ny,*)
     complex(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_C2C_cpu_cpp
   end function gather_fft_3D_C2C_cpu_cpp
   function gather_fft_3D_D2Z_cpu_cpp(nx,ny,nz,fxyz,fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     real(dp) :: fxyz(nx,ny,*)
     complex(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_D2Z_cpu_cpp
   end function gather_fft_3D_D2Z_cpu_cpp
   function gather_fft_3D_S2C_cpu_cpp(nx,ny,nz,s_fxyz,s_fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     real(sp) :: s_fxyz(nx,ny,*)
     complex(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_S2C_cpu_cpp
   end function gather_fft_3D_S2C_cpu_cpp
   function gather_fft_3D_Z2D_cpu_cpp(nx,ny,nz,fuvw,fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     complex(dp) :: fuvw(nx,ny,*)
     real(dp) :: fhpq(nx,ny,*)
     integer :: gather_fft_3D_Z2D_cpu_cpp
   end function gather_fft_3D_Z2D_cpu_cpp
   function gather_fft_3D_C2S_cpu_cpp(nx,ny,nz,s_fuvw,s_fhpq)
     use simple_defs
     integer  :: nx,ny,nz
     complex(sp) :: s_fuvw(nx,ny,*)
     real(sp) :: s_fhpq(nx,ny,*)
     integer :: gather_fft_3D_C2S_cpu_cpp
   end function gather_fft_3D_C2S_cpu_cpp

#endif

end interface

contains 
  !CONSTRUCTORS
  !> \brief is a fft123D constructor
  function constructor_fft_cpu(nx,ny,nz,fft_in) result(c_devFFT_cpu)
    type(fft_cpu) :: c_devFFT_cpu
    integer,intent(inout),optional :: nx,ny,nz
    integer,intent(in) :: fft_in

    if ( present(nx) ) then
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_cpu%new_fft_1D_DZ_cpu(nx)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_cpu%new_fft_1D_ZZ_cpu(nx)
       end if
    else if ( present(nx) .and. present(ny) ) then
       !single precision
       if ( fft_in == FFT_S2C .or. fft_in == FFT_C2S ) then
          call c_devFFT_cpu%new_fft_2D_SC_cpu(nx,ny)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_cpu%new_fft_2D_CC_cpu(nx,ny)
       end if
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_cpu%new_fft_2D_DZ_cpu(nx,ny)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_cpu%new_fft_2D_ZZ_cpu(nx,ny)
       end if
    else if ( present(nx) .and. present(ny) .and. present(nz) ) then
       !single precision
       if ( fft_in == FFT_S2C .or. fft_in == FFT_C2S ) then
          call c_devFFT_cpu%new_fft_3D_SC_cpu(nx,ny,nz)
       else if ( fft_in == FFT_C2C ) then
          call c_devFFT_cpu%new_fft_3D_CC_cpu(nx,ny,nz)
       end if
       !double precision
       if ( fft_in == FFT_D2Z .or. fft_in == FFT_Z2D ) then
          call c_devFFT_cpu%new_fft_3D_DZ_cpu(nx,ny,nz)
       else if ( fft_in == FFT_Z2Z ) then
          call c_devFFT_cpu%new_fft_3D_ZZ_cpu(nx,ny,nz)
       end if
    end if

  end function constructor_fft_cpu

  ! 3D constructor method 

  subroutine new_fft_3D_CC_cpu(c_devFFT_cpu,nx,ny,nz)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny,nz

    !complex(sp)
    allocate( c_devFFT_cpu%s_u(nx))
    allocate( c_devFFT_cpu%s_v(ny))
    allocate( c_devFFT_cpu%s_w(nz))
    allocate( c_devFFT_cpu%s_fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_3D_CC_cpu
  subroutine new_fft_3D_SC_cpu(c_devFFT_cpu,nx,ny,nz)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny,nz

    !real(sp)
    allocate( c_devFFT_cpu%s_x(nx))
    allocate( c_devFFT_cpu%s_y(ny))
    allocate( c_devFFT_cpu%s_z(nz))
    allocate(c_devFFT_cpu%fxyz(nx,ny,nz))
    !complex(sp)
    allocate( c_devFFT_cpu%s_u(nx))
    allocate( c_devFFT_cpu%s_v(ny))
    allocate( c_devFFT_cpu%s_w(nz))
    allocate(c_devFFT_cpu%s_fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_3D_SC_cpu

  subroutine new_fft_3D_ZZ_cpu(c_devFFT_cpu,nx,ny,nz)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny,nz

    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate( c_devFFT_cpu%v(ny))
    allocate( c_devFFT_cpu%w(nz))
    allocate(c_devFFT_cpu%fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_3D_ZZ_cpu
  subroutine new_fft_3D_DZ_cpu(c_devFFT_cpu,nx,ny,nz)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny,nz

    !real(dp)
    allocate( c_devFFT_cpu%x(nx))
    allocate( c_devFFT_cpu%y(ny))
    allocate( c_devFFT_cpu%z(nz))
    allocate(c_devFFT_cpu%fxyz(nx,ny,nz))
    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate( c_devFFT_cpu%v(ny))
    allocate( c_devFFT_cpu%w(nz))
    allocate(c_devFFT_cpu%fuvw(nx,ny,nz))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_3D_DZ_cpu
  ! 2D constructor method 

  subroutine new_fft_2D_CC_cpu(c_devFFT_cpu,nx,ny)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_cpu%s_u(nx))
    allocate( c_devFFT_cpu%s_v(ny))
    allocate( c_devFFT_cpu%s_fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_2D_CC_cpu
  subroutine new_fft_2D_SC_cpu(c_devFFT_cpu,nx,ny)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices

    !real(dp)
    allocate( c_devFFT_cpu%s_x(nx))
    allocate( c_devFFT_cpu%s_y(ny))
    allocate( c_devFFT_cpu%s_fxy(nx,ny))
    !complex(dp)
    allocate( c_devFFT_cpu%s_u(nx))
    allocate( c_devFFT_cpu%s_v(ny))
    allocate( c_devFFT_cpu%s_fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_2D_SC_cpu

  subroutine new_fft_2D_ZZ_cpu(c_devFFT_cpu,nx,ny)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate( c_devFFT_cpu%v(ny))
    allocate(c_devFFT_cpu%fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_2D_ZZ_cpu
  subroutine new_fft_2D_DZ_cpu(c_devFFT_cpu,nx,ny)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx,ny
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices
    
    !real(dp)
    allocate( c_devFFT_cpu%x(nx))
    allocate( c_devFFT_cpu%y(ny))
    allocate(c_devFFT_cpu%fxy(nx,ny))
    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate( c_devFFT_cpu%v(ny))
    allocate(c_devFFT_cpu%fuv(nx,ny))

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_2D_DZ_cpu

  ! 1D constructor method 
  subroutine new_fft_1D_ZZ_cpu(c_devFFT_cpu,nx)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices

    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate(c_devFFT_cpu%fu(nx))

    call gather_fft(c_devFFT_cpu,t_devFFT)

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_1D_ZZ_cpu
  subroutine new_fft_1D_DZ_cpu(c_devFFT_cpu,nx)
    class(fft_cpu) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: nx
    !kill the existing object before allocating a new one
    call c_devFFT_cpu%kill_fft_cpu

    !allocating arrays for details for n devices
    
    !real(dp)
    allocate( c_devFFT_cpu%x(nx))
    allocate(c_devFFT_cpu%fx(nx))
    !complex(dp)
    allocate( c_devFFT_cpu%u(nx))
    allocate(c_devFFT_cpu%fu(nx))

    call gather_fft(c_devFFT_cpu,t_devFFT)

    !set the existence of the object to true
    c_devFFT_cpu%existence_c_devFFT_cpu = .true.

    return
  end subroutine new_fft_1D_DZ_cpu

  !WORKERS METHODS
  subroutine gather_fft(c_devFFT_cpu,t_devFFT)
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    !TODO: general gather method to parse general info may be removed later
    return
  end subroutine gather_fft

  ! 3D data
  !getting the fourier transform on CPU for
  ! Z2Z
  subroutine gather_fft_3D_Z2Z_cpu(c_devFFT_cpu,nx,ny,nz,fuvw,fhpq,sign)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz,sign
    complex(dp) :: fuvw(nx,ny,*)
    complex(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_Z2Z_cpu_c
       function gather_fft_3D_Z2Z_cpu_c(nx,ny,nz,fuvw,fhpq,sign)
         use simple_defs
         integer  :: nx,ny,nz
         integer  :: sign
         complex(dp) :: fuvw(nx,ny,*)
         complex(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_Z2Z_cpu_c
       end function gather_fft_3D_Z2Z_cpu_c
    end interface external_c_function_gather_fft_3D_Z2Z_cpu_c
    
#if defined (MACOSX)
    rc = gather_fft_3D_Z2Z_cpu_c(nx,ny,nz,fuvw,fhpq,sign)
#elif defined (LINUX)
    rc = gather_fft_3D_Z2Z_cpu_cpp(nx,ny,nz,fuvw,fhpq,sign)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_3D_Z2Z_cpu
  ! C2C
  subroutine gather_fft_3D_C2C_cpu(c_devFFT_cpu,nx,ny,nz,s_fuvw,s_fhpq,sign)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz,sign
    complex(sp) :: s_fuvw(nx,ny,*)
    complex(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_C2C_cpu_c
       function gather_fft_3D_C2C_cpu_c(nx,ny,nz,s_fuvw,s_fhpq,sign)
         use simple_defs
         integer  :: nx,ny,nz
         integer  :: sign
         complex(sp) :: s_fuvw(nx,ny,*)
         complex(sp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_C2C_cpu_c
       end function gather_fft_3D_C2C_cpu_c
    end interface external_c_function_gather_fft_3D_C2C_cpu_c

#if defined (MACOSX)
    rc = gather_fft_3D_C2C_cpu_c(nx,ny,nz,s_fuvw,s_fhpq,sign)
#elif defined (LINUX)
    rc = gather_fft_3D_C2C_cpu_cpp(nx,ny,nz,s_fuvw,s_fhpq,sign)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_3D_C2C_cpu
  ! D2Z
  subroutine gather_fft_3D_D2Z_cpu(c_devFFT_cpu,nx,ny,nz,fxyz,fhpq)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz
    real(dp) :: fxyz(nx,ny,*)
    complex(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_D2Z_cpu_c
       function gather_fft_3D_D2Z_cpu_c(nx,ny,nz,fxyz,fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         real(dp) :: fxyz(nx,ny,*)
         complex(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_D2Z_cpu_c
       end function gather_fft_3D_D2Z_cpu_c
    end interface external_c_function_gather_fft_3D_D2Z_cpu_c
#if defined (MACOSX)
    rc = gather_fft_3D_D2Z_cpu_c(nx,ny,nz,fxyz,fhpq)
#elif defined (LINUX)
    rc = gather_fft_3D_D2Z_cpu_cpp(nx,ny,nz,fxyz,fhpq)
#else
    call c_devFFT_cpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_D2Z_cpu
  ! S2C
  subroutine gather_fft_3D_S2C_cpu(c_devFFT_cpu,nx,ny,nz,s_fxyz,s_fhpq)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz
    real(sp) :: s_fxyz(nx,ny,*)
    complex(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_S2C_cpu_c
       function gather_fft_3D_S2C_cpu_c(nx,ny,nz,s_fxyz,s_fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         real(sp) :: s_fxyz(nx,ny,*)
         complex(sp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_S2C_cpu_c
       end function gather_fft_3D_S2C_cpu_c
    end interface external_c_function_gather_fft_3D_S2C_cpu_c
#if defined (MACOSX)
    rc = gather_fft_3D_S2C_cpu_c(nx,ny,nz,s_fxyz,s_fhpq)
#elif defined (LINUX)
    rc = gather_fft_3D_S2C_cpu_cpp(nx,ny,nz,s_fxyz,s_fhpq)
#else
    call c_devFFT_cpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_S2C_cpu
  ! Z2D
  subroutine gather_fft_3D_Z2D_cpu(c_devFFT_cpu,nx,ny,nz,fuvw,fhpq)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz
    complex(dp) :: fuvw(nx,ny,*)
    real(dp) :: fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_Z2D_cpu_c
       function gather_fft_3D_Z2D_cpu_c(nx,ny,nz,fuvw,fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         complex(dp) :: fuvw(nx,ny,*)
         real(dp) :: fhpq(nx,ny,*)
         integer :: gather_fft_3D_Z2D_cpu_c
       end function gather_fft_3D_Z2D_cpu_c
    end interface external_c_function_gather_fft_3D_Z2D_cpu_c
#if defined (MACOSX)
    rc = gather_fft_3D_Z2D_cpu_c(nx,ny,nz,fuvw,fhpq)
#elif defined (LINUX)
    rc = gather_fft_3D_Z2D_cpu_cpp(nx,ny,nz,fuvw,fhpq)
#else
    call c_devFFT_cpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_Z2D_cpu
  ! C2S
  subroutine gather_fft_3D_C2S_cpu(c_devFFT_cpu,nx,ny,nz,s_fuvw,s_fhpq)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,nz
    complex(sp) :: s_fuvw(nx,ny,*)
    real(sp) :: s_fhpq(nx,ny,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_3D_C2S_cpu_c
       function gather_fft_3D_C2S_cpu_c(nx,ny,nz,s_fuvw,s_fhpq)
         use simple_defs
         integer  :: nx,ny,nz
         complex(sp) :: s_fuvw(nx,ny,*)
         real(sp) :: s_fhpq(nx,ny,*)
         integer :: gather_fft_3D_C2S_cpu_c
       end function gather_fft_3D_C2S_cpu_c
    end interface external_c_function_gather_fft_3D_C2S_cpu_c
#if defined (MACOSX)
    rc = gather_fft_3D_C2S_cpu_c(nx,ny,nz,s_fuvw,s_fhpq)
#elif defined (LINUX)
    rc = gather_fft_3D_C2S_cpu_cpp(nx,ny,nz,s_fuvw,s_fhpq)
#else
    call c_devFFT_cpu%terminate()
#endif
    return
  end subroutine gather_fft_3D_C2S_cpu

  ! 2D data
  !getting the fourier transform on CPU for
  ! Z2Z
  subroutine gather_fft_2D_Z2Z_cpu(c_devFFT_cpu,nx,ny,fuv,fhp,sign)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,sign
    complex(dp) :: fuv(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_Z2Z_cpu_c
       function gather_fft_2D_Z2Z_cpu_c(nx, ny,fuv,fhp,sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         complex(dp) :: fuv(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_Z2Z_cpu_c
       end function gather_fft_2D_Z2Z_cpu_c
    end interface external_c_function_gather_fft_2D_Z2Z_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_Z2Z_cpu_c(nx, ny, fuv, fhp, sign)
#elif defined (LINUX)
    rc = gather_fft_2D_Z2Z_cpu_cpp(nx, ny,fuv,fhp,sign)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_Z2Z_cpu
  ! C2C
  subroutine gather_fft_2D_C2C_cpu(c_devFFT_cpu,nx,ny,s_fuv,s_fhp,sign)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny,sign
    complex(sp) :: s_fuv(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_C2C_cpu_c
       function gather_fft_2D_C2C_cpu_c(nx, ny,s_fuv,s_fhp,sign)
         use simple_defs
         integer  :: nx,ny
         integer  :: sign
         complex(sp) :: s_fuv(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_C2C_cpu_c
       end function gather_fft_2D_C2C_cpu_c
    end interface external_c_function_gather_fft_2D_C2C_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_C2C_cpu_c(nx, ny, s_fuv, s_fhp, sign)
#elif defined (LINUX)
    rc = gather_fft_2D_C2C_cpu_cpp(nx, ny,s_fuv,s_fhp,sign)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_C2C_cpu
  ! D2Z
  subroutine gather_fft_2D_D2Z_cpu(c_devFFT_cpu,nx,ny,fxy,fhp)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny
    real(dp) :: fxy(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_D2Z_cpu_c
       function gather_fft_2D_D2Z_cpu_c(nx, ny,fxy,fhp)
         use simple_defs
         integer  :: nx,ny
         real(dp) :: fxy(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_D2Z_cpu_c
       end function gather_fft_2D_D2Z_cpu_c
    end interface external_c_function_gather_fft_2D_D2Z_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_D2Z_cpu_c(nx, ny, fxy, fhp)
#elif defined (LINUX)
    rc = gather_fft_2D_D2Z_cpu_cpp(nx, ny,fxy,fhp)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_D2Z_cpu
  ! S2C
  subroutine gather_fft_2D_S2C_cpu(c_devFFT_cpu,nx,ny,s_fxy,s_fhp)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny
    real(sp) :: s_fxy(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_S2C_cpu_c
       function gather_fft_2D_S2C_cpu_c(nx, ny,s_fxy,s_fhp)
         use simple_defs
         integer  :: nx,ny
         real(sp) :: s_fxy(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_S2C_cpu_c
       end function gather_fft_2D_S2C_cpu_c
    end interface external_c_function_gather_fft_2D_S2C_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_S2C_cpu_c(nx, ny, s_fxy, s_fhp)
#elif defined (LINUX)
    rc = gather_fft_2D_S2C_cpu_cpp(nx, ny,s_fxy,s_fhp)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_S2C_cpu
  ! Z2D
  subroutine gather_fft_2D_Z2D_cpu(c_devFFT_cpu,nx,ny,fhp,fxy)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny
    real(dp) :: fxy(nx,*)
    complex(dp) :: fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_Z2D_cpu_c
       function gather_fft_2D_Z2D_cpu_c(nx, ny,fhp,fxy)
         use simple_defs
         integer  :: nx,ny
         real(dp) :: fxy(nx,*)
         complex(dp) :: fhp(nx,*)
         integer :: gather_fft_2D_Z2D_cpu_c
       end function gather_fft_2D_Z2D_cpu_c
    end interface external_c_function_gather_fft_2D_Z2D_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_Z2D_cpu_c(nx, ny, fhp, fxy)
#elif defined (LINUX)
    rc = gather_fft_2D_Z2D_cpu_cpp(nx, ny,fhp,fxy)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_Z2D_cpu
  ! C2S
  subroutine gather_fft_2D_C2S_cpu(c_devFFT_cpu,nx,ny,s_fhp,s_fxy)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,ny
    real(sp) :: s_fxy(nx,*)
    complex(sp) :: s_fhp(nx,*)
    integer :: rc !return code
    interface external_c_function_gather_fft_2D_C2S_cpu_c
       function gather_fft_2D_C2S_cpu_c(nx, ny,s_fhp,s_fxy)
         use simple_defs
         integer  :: nx,ny
         real(sp) :: s_fxy(nx,*)
         complex(sp) :: s_fhp(nx,*)
         integer :: gather_fft_2D_C2S_cpu_c
       end function gather_fft_2D_C2S_cpu_c
    end interface external_c_function_gather_fft_2D_C2S_cpu_c

#if defined (MACOSX)
    rc = gather_fft_2D_C2S_cpu_c(nx, ny, s_fhp, s_fxy)
#elif defined (LINUX)
    rc = gather_fft_2D_C2S_cpu_cpp(nx, ny,s_fhp,s_fxy)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_2D_C2S_cpu

  ! 1D data
  !getting the fourier transform on CPU for
  ! Z2Z
  subroutine gather_fft_1D_Z2Z_cpu(c_devFFT_cpu,nx,fu,fh,sign)
    use simple_defs
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer  :: nx,sign
    complex(dp) :: fu(*)
    complex(dp) :: fh(*)
    integer :: rc !return code
    interface external_c_function_gather_fft_1D_Z2Z_cpu_c
       function gather_fft_1D_Z2Z_cpu_c(nx,fu,fh,sign)
         use simple_defs
         integer  :: nx
         integer  :: sign
         complex(dp) :: fu(*)
         complex(dp) :: fh(*)
         integer :: gather_fft_1D_Z2Z_cpu_c
       end function gather_fft_1D_Z2Z_cpu_c
    end interface external_c_function_gather_fft_1D_Z2Z_cpu_c

#if defined (MACOSX)
    rc = gather_fft_1D_Z2Z_cpu_c(nx,fu,fh,sign)
#elif defined (LINUX)
    rc = gather_fft_1D_Z2Z_cpu_cpp(nx,fu,fh,sign)
#else
    call c_devFFT_cpu%terminate()
#endif

    return
  end subroutine gather_fft_1D_Z2Z_cpu

  !SETTERS
  ! 1D
  subroutine set_Ddata_1D(c_devFFT_cpu,nx,x,data1D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer  :: nx
    real(dp) :: x(*)
    real(dp) :: data1D(*)

    c_devFFT_cpu%x(1:nx) = x(1:nx)
    c_devFFT_cpu%fx(1:nx) = data1D(1:nx)

    return
  end subroutine set_Ddata_1D
  subroutine set_Zdata_1D(c_devFFT_cpu,nx,u,Zdata1D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer  :: nx
    complex(dp) :: u(*)
    complex(dp) :: Zdata1D(*)

    c_devFFT_cpu%u(1:nx) = u(1:nx)
    c_devFFT_cpu%fu(1:nx) = Zdata1D(1:nx)

    return
  end subroutine set_Zdata_1D
  ! 2D
  subroutine set_Ddata_2D(c_devFFT_cpu,lda,data2D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer :: lda
    real(dp) :: data2D(lda,*)
    return
  end subroutine set_Ddata_2D
  subroutine set_Zdata_2D(c_devFFT_cpu,nx,ny,Zdata2D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer :: nx,ny
    complex(dp) :: Zdata2D(nx,*)
    c_devFFT_cpu%fuv(1:nx,1:ny) = Zdata2D(1:nx,1:ny)
    return
  end subroutine set_Zdata_2D
  ! 3D
  subroutine set_Ddata_3D(c_devFFT_cpu,lda,data3D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer :: lda
    real(dp) :: data3D(lda,*)
    return
  end subroutine set_Ddata_3D
  subroutine set_Zdata_3D(c_devFFT_cpu,nx,ny,nz,Zdata3D)
    implicit none
    class(fft_cpu),intent(inout) :: c_devFFT_cpu
    integer :: nx,ny,nz
    complex(dp) :: Zdata3D(nx,ny,*)
    c_devFFT_cpu%fuvw(1:nx,1:ny,1:nz) = Zdata3D(1:nx,1:ny,1:nz)
    return
  end subroutine set_Zdata_3D

  !GETTERS
  !double precision
  function get_Ddata_1D(c_devFFT_cpu,nx) result(Ddata1D_out)
    class(fft_cpu),intent(in) :: c_devFFT_cpu
    integer :: nx
    real(dp) :: Ddata1D_out(1:nx)
    Ddata1D_out(1:nx) = c_devFFT_cpu%fx(1:nx)
  end function get_Ddata_1D
  function get_x_1D(c_devFFT_cpu,nx) result(x_out)
    class(fft_cpu),intent(in) :: c_devFFT_cpu
    integer :: nx
    real(dp) :: x_out(1:nx)
    x_out(1:nx) = c_devFFT_cpu%x(1:nx)
  end function get_x_1D
  ! complex double precision
  function get_Zdata_1D(c_devFFT_cpu,nx) result(Zdata1D_out)
    class(fft_cpu),intent(in) :: c_devFFT_cpu
    integer :: nx
    complex(dp) :: Zdata1D_out(1:nx)
    Zdata1D_out(1:nx) = c_devFFT_cpu%fu(1:nx)
  end function get_Zdata_1D
  function get_u_1D(c_devFFT_cpu,nx) result(u_out)
    class(fft_cpu),intent(in) :: c_devFFT_cpu
    integer :: nx
    complex(dp) :: u_out(1:nx)
    u_out(1:nx) = c_devFFT_cpu%u(1:nx)
  end function get_u_1D

  !INITIALISERS
  subroutine initialise_fft(c_devFFT_cpu,idev,t_devFFT)
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    type(fft_devType) :: t_devFFT
    integer :: idev
    integer :: rc !return code
    interface external_c_function_init
    end interface external_c_function_init

#if defined (MACOSX)
#elif defined (LINUX)
#else
    call c_devFFT_cpu%terminate()
#endif
    return
  end subroutine initialise_fft

  !GREETERS
  !> \brief hello greeting routine for the object
  subroutine hello_fft_cpu(err)
    implicit none
    integer :: err
    !start of the execution commands
    write(*,*) "Hello fft CPU world"
    write(*,*)
    return
  end subroutine hello_fft_cpu

  !> \brief bye greeting routine for the object
  subroutine bye_fft_cpu()
    implicit none
    !start of the execution commands
    write(*,*) "Bye fft CPU world"
    write(*,*)
    return
  end subroutine bye_fft_cpu

  !TERMINATORS because of incompatiblilties errors
  subroutine terminate(c_devFFT_cpu)
    class(fft_cpu) :: c_devFFT_cpu
    call c_devFFT_cpu%kill_fft_cpu
    write(*,*)"**************************WARNING*******************************"
    write(*,*)"*You are trying to run the the app from *outer space* other    *"
    write(*,*)"*where the dimensions and the real physics operates!!!         *"
    write(*,*)"*Maybe trying to compile on a Linux or MacOSX operating system *"
    write(*,*)"*will bing you back to reality and get you some nice results   *"
    write(*,*)"****************************************************************"
    call bye_fft_cpu()
    stop
    return
  end subroutine terminate

  !DESTRUCTOR
  !> \brief is a fft_cpu destructor
  subroutine kill_fft_cpu(c_devFFT_cpu)
    class(fft_cpu), intent(inout) :: c_devFFT_cpu
    if ( c_devFFT_cpu%existence_c_devFFT_cpu) then

       if(allocated(c_devFFT_cpu%x )) deallocate(c_devFFT_cpu%x )
       if(allocated(c_devFFT_cpu%y )) deallocate(c_devFFT_cpu%y )
       if(allocated(c_devFFT_cpu%z )) deallocate(c_devFFT_cpu%z )

       if(allocated(c_devFFT_cpu%fx)) deallocate(c_devFFT_cpu%fx)
       if(allocated(c_devFFT_cpu%fxy)) deallocate(c_devFFT_cpu%fxy)
       if(allocated(c_devFFT_cpu%fxyz)) deallocate(c_devFFT_cpu%fxyz)

       if(allocated(c_devFFT_cpu%u )) deallocate(c_devFFT_cpu%u )
       if(allocated(c_devFFT_cpu%v )) deallocate(c_devFFT_cpu%v )
       if(allocated(c_devFFT_cpu%w )) deallocate(c_devFFT_cpu%w )

       if(allocated(c_devFFT_cpu%fu)) deallocate(c_devFFT_cpu%fu)
       if(allocated(c_devFFT_cpu%fuv)) deallocate(c_devFFT_cpu%fuv)
       if(allocated(c_devFFT_cpu%fuvw)) deallocate(c_devFFT_cpu%fuvw)

       c_devFFT_cpu%existence_c_devFFT_cpu = .false.

    end if
    return
  end subroutine kill_fft_cpu

end module fft123D_cpu
