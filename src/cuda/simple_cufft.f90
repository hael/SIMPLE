!>  \brief  CUFFT NVIDIA interfaces
! based on cufft_m.cuf CUDA_Fortran examples
module simple_cufft

  integer, parameter, public :: CUFFT_FORWARD = -1
  integer, parameter, public :: CUFFT_INVERSE = 1
  integer, parameter, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
  integer, parameter, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
  integer, parameter, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
  integer, parameter, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
  integer, parameter, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
  integer, parameter, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex
  ! ------------
  ! cufftDestroy
  ! ------------

  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftDestroy
#endif
       use iso_c_binding
       integer(c_int),value:: plan
     end subroutine cufftDestroy
  end interface cufftDestroy

  ! A generic interface such as "cufftExec" which contained all the R2C, C2C, ...
  ! subroutine variants would work fine as long as transforms were not done in place.
  ! The types and kinds associated with idata and odata would determine the
  ! correct routine.
  !
  ! However, this appoach would break down when performing in-place transforms.
  ! For example, suppose one wanted to perform a real to complex transform in-place.
  ! The soubroutine call:
  !    call cufftExec(plan, data, data)
  ! would result in a compile time error if "data" were of type real or
  ! the cufftExecC2C would be executed if "data" were of type complex

  ! So, we define a separate interface for each R2C, C2C, ... variant and
  ! use the "!dir$ ignore_tkr" directive to ignore compile-time checks of
  ! data type, kind, and ranks.

  ! ----------------------------------
  ! cufftExec[C2C|R2C|C2R|Z2Z|Z2D|D2Z]
  ! ----------------------------------

  interface cufftExecC2C
     subroutine cufftExecC2C(plan, idata, odata, direction) &
          & bind(C,name='cufftExecC2C')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecC2C
#endif
       use iso_c_binding
       use simple_defs
       integer(c_int), value:: plan
       !$ ignore_tkr idata, odata
       complex(sp), device :: idata(*), odata(*)
       integer(c_int), value:: direction
     end subroutine cufftExecC2C
  end interface cufftExecC2C

  interface cufftExecR2C
     subroutine cufftExecR2C(plan, idata, odata) bind(C,name='cufftExecR2C')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecR2C
#endif
       use iso_c_binding
       use simple_defs
       integer (c_int), value :: plan
       !$ ignore_tkr idata, odata
       real (sp), device :: idata(*)
       complex (sp), device :: odata(*)
     end subroutine cufftExecR2C
  end interface cufftExecR2C

  interface cufftExecC2R
     subroutine cufftExecC2R(plan, idata, odata) bind(C,name='cufftExecC2R')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecC2R
#endif
       use iso_c_binding
       use simple_defs
       integer (c_int), value :: plan
       !$ ignore_tkr idata, odata
       complex (sp), device :: idata(*)
       real (sp), device :: odata(*)
     end subroutine cufftExecC2R
  end interface cufftExecC2R

  interface cufftExecZ2Z
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          & bind(C,name='cufftExecZ2Z')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecZ2Z
#endif
       use iso_c_binding
       use simple_defs
       integer(c_int),value:: plan
       !$ ignore_tkr idata, odata
       complex(dp), device:: idata(*), odata(*)
       integer(c_int),value:: direction
     end subroutine cufftExecZ2Z
  end interface cufftExecZ2Z

  interface cufftExecD2Z
     subroutine cufftExecD2Z(plan, idata, odata) bind(C,name='cufftExecD2Z')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecD2Z
#endif
       use iso_c_binding
       use simple_defs
       integer (c_int), value :: plan
       !$ ignore_tkr idata, odata
       real (dp), device :: idata(*)
       complex (dp), device :: odata(*)
     end subroutine cufftExecD2Z
  end interface cufftExecD2Z

  interface cufftExecZ2D
     subroutine cufftExecZ2D(plan, idata, odata) bind(C,name='cufftExecZ2D')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftExecZ2D
#endif
       use iso_c_binding
       use simple_defs
       integer (c_int), value :: plan
       !$ ignore_tkr idata, odata
       complex (dp), device :: idata(*)
       real (dp), device :: odata(*)
     end subroutine cufftExecZ2D
  end interface cufftExecZ2D

  interface cufftExec
     subroutine cufftExec(plan, transform, idata, odata, direction)
       integer :: plan, transform
       !$ ignore_tkr idata, odata
       real, device :: idata(*), odata(*)
       integer, optional :: direction
     end subroutine cufftExec
  end interface cufftExec

  ! -----------------
  ! cufftPlan[1|2|3]d
  ! -----------------

  interface cufftPlan1d
     subroutine cufftPlan1d(plan, nx, type, batch) bind(C,name='cufftPlan1d')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftPlan1d
#endif
       use iso_c_binding
       integer(c_int):: plan
       integer(c_int),value:: nx, batch,type
     end subroutine cufftPlan1d
  end interface cufftPlan1d

  interface cufftPlan2d
     subroutine cufftPlan2d(plan, nx, ny, type) bind(C,name='cufftPlan2d')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftPlan2d
#endif
       use iso_c_binding
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, type
     end subroutine cufftPlan2d
  end interface cufftPlan2d

  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type) bind(C,name='cufftPlan3d')
#ifdef _WIN32
       !dec$ attributes stdcall, decorate :: cufftPlan3d
#endif
       use iso_c_binding
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, nz, type
     end subroutine cufftPlan3d
  end interface cufftPlan3d

  ! interface cufftPlanMany
  !    subroutine cufftPlanMany(plan, rank, n, inembed, istride, idist, &
  !         onembed, ostride, odist,  &
  !         type, batch) bind(C,name='cufftPlanMany')
  !      use iso_c_binding
  !      implicit none
  !      !$ ignore_tkr n, inembed, onembed
  !      type(c_ptr) :: plan
  !      integer(c_int) :: n, inembed, onembed
  !      integer(c_int), value:: rank, istride, ostride, idist, odist, type, batch
  !    end subroutine cufftPlanMany
  ! end interface cufftPlanMany

  ! interface cufftPlan2d
  !    module procedure cufftPlan2Dswap
  ! end interface cufftPlan2d

  ! interface cufftPlan2dC
  !    subroutine cufftPlan2d(plann, nxn, nyn, typen) &
  !         bind(C,name='cufftPlan2d')
  !      use iso_c_binding
  !      type(c_ptr):: plann
  !      integer(c_int),value:: nxn, nyn, typen
  !    end subroutine cufftPlan2d
  ! end interface cufftPlan2dC

! contains

!   subroutine cufftPlan2Dswap(plan,nx,ny, type)
!     use iso_c_binding
!     type(c_ptr):: plan
!     integer(c_int),value:: nx, ny, type
!     call cufftPlan2dC(plan,ny,nx,type)
!   end subroutine cufftPlan2Dswap

end module simple_cufft

! A general subroutine that uses the transform type passed in "transform" to
! determine which routine to call.  The interfaces for the routines specified
! above use the "!dir$ ignore_tkr" directive to remove checking of
! type, kind, and rank to facilitate certain in-place transforms (eg. cufftExecR2C)
! and multidimensional arrays

subroutine cufftExec(plan, transform, idata, odata, direction)
  use simple_cufft
  implicit none

  integer :: plan, transform
  real, device :: idata(*), odata(*)
  integer, optional :: direction


  select case (transform)
  case (CUFFT_R2C)
     call cufftExecR2C(plan, idata, odata)
  case (CUFFT_C2R)
     call cufftExecC2R(plan, idata, odata)
  case (CUFFT_C2C)
     if (.not. present(direction)) then
        write(*,*) 'Error in cufft: C2C transform called without specifying direction'
        stop
     end if
     call cufftExecC2C(plan, idata, odata, direction)
  case (CUFFT_D2Z)
     call cufftExecD2Z(plan, idata, odata)
  case (CUFFT_Z2D)
     call cufftExecZ2D(plan, idata, odata)
  case (CUFFT_Z2Z)
     if (.not. present(direction)) then
        write(*,*) 'Error in cufft: Z2Z transform called without specifying direction'
        stop
     end if
     call cufftExecZ2Z(plan, idata, odata, direction)
  case default
     write(*,*) 'Invalid transform type passed to cufftExec'
  end select
end subroutine cufftExec


! FFTW function  |                                  CUFFT function
!-----------------------------------------------------------------
!   fftw_plan_dft_1d(),
!     fftw_plan_dft_r2c_1d(),                      cufftPlan1d()
!     fftw_plan_dft_c2r_1d()
!-----------------------------------------------------------------
!   fftw_plan_dft_2d(),
!     fftw_plan_dft_r2c_2d(),
!     fftw_plan_dft_c2r_2d()                       cufftPlan2d()
!-----------------------------------------------------------------
!   fftw_plan_dft_3d(),
!   fftw_plan_dft_r2c_3d(),
!   fftw_plan_dft_c2r_3d()                         cufftPlan3d()
!-----------------------------------------------------------------
!   fftw_plan_dft(),        \
!   fftw_plan_dft_r2c(),     |
!   fftw_plan_dft_c2r()     /                       cufftPlanMany()
!-----------------------------------------------------------------
!   fftw_plan_many_dft(),     \
!   fftw_plan_many_dft_r2c(),  |                   cufftPlanMany()
!   fftw_plan_many_dft_c2r()  /
!-----------------------------------------------------------------
!   fftw_execute()                                 cufftExecC2C(), cufftExecZ2Z(),
!                                                  cufftExecR2C(), cufftExecD2Z(),
!                                                  cufftExecC2R(), cufftExecZ2D()
!-----------------------------------------------------------------
!   fftw_destroy_plan()                            cufftDestroy()
