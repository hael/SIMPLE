
!===============================================================================
! Copyright 2011-2017 Intel Corporation All Rights Reserved.
!
! The source code,  information  and material  ("Material") contained  herein is
! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
! Material  contains  proprietary  information  of  Intel or  its suppliers  and
! licensors.  The Material is protected by  worldwide copyright  laws and treaty
! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
! in any way without Intel's prior express written permission.  No license under
! any patent,  copyright or other  intellectual property rights  in the Material
! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
! property rights must be express and approved by Intel in writing.
!
! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
! suppliers or licensors in any way.
!===============================================================================

! Content:
! An example of using Intel(R) MKL DFTI configuration parameter
! DFTI_CONJUGATE_EVEN_STORAGE. The parameter defines layout of complex data in
! the backward domain of real-to-complex FFT.
!
! Values:
! DFTI_COMPLEX_REAL (default for 1d and 2d transforms, not recommended)
!     represent the complex data by real and imaginary parts packed in a real
!     array as defined by DFTI_PACKED_FORMAT configuration parameter.
!     This example shows usage of this default setting.
!
! DFTI_COMPLEX_COMPLEX (recommented, default for 3d and higher rank FFTs)
!     represent the complex data by complex elements. For the recommended use
!     of DFTI_CONJUGATE_EVEN_STORAGE see examples of real DFT.
!
! Note: DFTI_COMPLEX_COMPLEX is recommended value for
!       DFTI_CONJUGATE_EVEN_STORAGE configuration.
!
!*****************************************************************************
module config_conjugate_even_storage
include 'simple_lib.f08'
implicit none
public :: test_CCE

integer, parameter :: WP = selected_real_kind(15,307)

contains

    subroutine test_CCE
#ifdef INTEL
        use MKL_DFTI, forget => DFTI_DOUBLE, DFTI_DOUBLE => DFTI_DOUBLE_R

  ! Sizes of 2D transform
  integer, parameter :: N1 = 8
  integer, parameter :: N2 = 14

  ! Arbitrary harmonic used to verify FFT
  integer, parameter :: H1 = 1
  integer, parameter :: H2 = -1

  ! Need double precision
  integer, parameter :: WP = selected_real_kind(15,307)

  ! Execution status
  integer :: status = 0, ignored_status

  ! Data array for in-place real FFT
  real(WP), allocatable:: x(:,:)

  ! DFTI descriptor handle
  type(DFTI_DESCRIPTOR), POINTER :: hand

  ! Strides define the data layout for forward and backward domains
  integer :: strides(3)



  hand => null()

  print *,"Example config_conjugate_even_storage"
  print *,"Real-to-complex in-place 2D FFT"
  print *,"Configuration parameters:"
  print *,"DFTI_PRECISION              = DFTI_DOUBLE"
  print *,"DFTI_FORWARD_DOMAIN         = DFTI_REAL"
  print *,"DFTI_DIMENSION              = 2"
  print '(" DFTI_LENGTHS                = /"I0","I0"/" )', N1, N2
  print *,"DFTI_CONJUGATE_EVEN_STORAGE = DFTI_COMPLEX_REAL"


  print *,"Create DFTI descriptor for 2D real in-place FFT"
  status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_REAL, 2, [N1,N2])
  if (0 /= status) goto 999

  ! This may be skipped because COMPLEX_REAL is default for 2D real FFT
  print *,"Set conjugate-even storage layout"
  status = DftiSetValue(hand, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL)
  if (0 /= status) goto 999

  print *,"Allocate array x(N1+2,N2+2), will suffice for all packed formats"
  allocate ( x(N1+2,N2+2), STAT = status)
  if (0 /= status) goto 999

  strides = [ 0, 1, size(x,dim=1) ]

  print '(" Set input  strides = ["3(I0:", ")"]")', strides
  status = DftiSetValue(hand, DFTI_INPUT_STRIDES, strides)
  if (0 /= status) goto 999

  print '(" Set output strides = ["3(I0:", ")"]")', strides
  status = DftiSetValue(hand, DFTI_OUTPUT_STRIDES, strides)
  if (0 /= status) goto 999


  print *,"===== Configure descriptor for CCS format ====="

  ! This may be skipped because CCS is default for 2D real FFT
  print *,"Set DFTI_PACKED_FORMAT = DFTI_CCS_FORMAT"
  status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_CCS_FORMAT)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize input data"
  call init(x, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result in CCS format"
  status = verificate(ccs, x, N1, N2, H1, H2)
  if (0 /= status) goto 999


  print *,"===== Configure descriptor for PACK format ====="

  print *,"Set DFTI_PACKED_FORMAT = DFTI_PACK_FORMAT"
  status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize input data"
  call init(x, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result in PACK format"
  status = verificate(pack, x, N1, N2, H1, H2)
  if (0 /= status) goto 999


  print *,"===== Configure descriptor for PERM format ====="

  print *,"Set DFTI_PACKED_FORMAT = DFTI_PERM_FORMAT"
  status = DftiSetValue(hand, DFTI_PACKED_FORMAT, DFTI_PERM_FORMAT)
  if (0 /= status) goto 999

  print *,"Commit DFTI descriptor"
  status = DftiCommitDescriptor(hand)
  if (0 /= status) goto 999

  print *,"Initialize input data"
  call init(x, N1, N2, H1, H2)

  print *,"Compute forward transform"
  status = DftiComputeForward(hand, x(:,1))
  if (0 /= status) goto 999

  print *,"Verify the result in PERM format"
  status = verificate(perm, x, N1, N2, H1, H2)
  if (0 /= status) goto 999

100 continue

  print *,"Release the DFTI descriptor"
  ignored_status = DftiFreeDescriptor(hand)

  if (allocated(x)) then
      print *,"Deallocate data arrays"
      deallocate(x)
  endif

  if (status == 0) then
    print *, "TEST PASSED"
    call exit(0)
  else
    print *, "TEST FAILED"
    call exit(1)
  end if

999 print '("  Error, status = ",I0)', status
  goto 100


#else
  print *,"Example config_conjugate_even_storage => Compiler /= INTEL"
  print *,"SKIPPING"
#endif

contains

  ! Compute mod(K*L,M) accurately
  pure real(WP) function moda(k,l,m)
    integer, intent(in) :: k,l,m
    integer*8 :: k8
    k8 = k
    moda = real(mod(k8*l,m),WP)
  end function moda

  ! Initialize x(:,:) to harmonic H
  subroutine init(x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)

    integer k1, k2
    real(WP), parameter:: TWOPI = 6.2831853071795864769_WP
    real(WP) :: factor

    factor = 2
    if (mod(2*(N1-H1),N1)==0 .and. mod(2*(N2-H2),N2)==0) factor = 1

    forall (k1=1:N1, k2=1:N2)
      x(k1,k2) = factor*cos(TWOPI*(moda(H1,k1-1,N1)/N1   &
    &    +                          moda(H2,k2-1,N2)/N2)) / (N1*N2)
    end forall
  end subroutine init

  ! Verify that x(:,:) has unit peak at (H1,H2)
  integer function verificate(unpack, x, N1, N2, H1, H2)
    integer N1, N2, H1, H2
    real(WP) :: x(:,:)
    external :: unpack

    integer k1, k2
    real(WP) err, errthr, maxerr
    complex(WP) res_exp, res_got

    ! Note, this simple error bound doesn't take into account error of
    ! input data
    errthr = 2.5 * log(real(N1*N2,WP)) / log(2.0_WP) * EPSILON(1.0_WP)
    print '("  Check if err is below errthr " G10.3)', errthr

    maxerr = 0.0_WP
    do k2 = 1, N2
      do k1 = 1, N1/2+1
        if (mod(k1-1-H1,N1)==0 .and. mod(k2-1-H2,N2)==0) then
          res_exp = 1.0_WP
        else if (mod(-k1+1-H1,N1)==0 .and. mod(-k2+1-H2,N2)==0) then
          res_exp = 1.0_WP
        else
          res_exp = 0.0_WP
        end if

        call unpack(res_got, x, size(x,dim=1), N1, N2, k1, k2)

        err = abs(res_got - res_exp)
        maxerr = max(err,maxerr)
        if (.not.(err < errthr)) then
          print '("  x("I0","I0"):"$)', k1,k2
          print '(" expected ("G24.17","G24.17"),"$)', res_exp
          print '(" got ("G24.17","G24.17"),"$)', res_got
          print '(" err "G10.3)', err
          print *,"  Verification FAILED"
          verificate = 100
          return
        end if
      end do
    end do
    print '("  Verified,  maximum error was " G10.3)', maxerr
    verificate = 0
  end function verificate

end subroutine test_CCE


! Fetch x(k1,k2) from the result of N1-by-N2 real FFT CCS-packed in x(:,:)
! Note: x should be embedded in a matrix at least x(N1+2,N2+2)
! Assume k1=1:N1, k2=1:N2
subroutine ccs(res, x, LD1, N1,N2, k1,k2)

  complex(WP) :: res
  integer :: LD1, N1, N2, k1, k2
  real(WP) :: x(LD1,*)

  real(WP) :: re, im

  if (k2 == 1) then
    if (k1 <= N1/2+1) then
      re =  x(1+2*(k1-1)+0, 1)
      im =  x(1+2*(k1-1)+1, 1)
    else
      re =  x(1+2*(N1-k1+1)+0, 1)
      im = -x(1+2*(N1-k1+1)+1, 1)
    end if
  else if (k1 == 1) then
    if (k2 <= N2/2+1) then
      re =  x(1, 1+2*(k2-1)+0)
      im =  x(1, 1+2*(k2-1)+1)
    else
      re =  x(1, 1+2*(N2-k2+1)+0)
      im = -x(1, 1+2*(N2-k2+1)+1)
    end if
  else if (k1-1 == N1-k1+1) then
    if (k2 <= N2/2+1) then
      re =  x(N1+1, 1+2*(k2-1)+0)
      im =  x(N1+1, 1+2*(k2-1)+1)
    else
      re =  x(N1+1, 1+2*(N2-k2+1)+0)
      im = -x(N1+1, 1+2*(N2-k2+1)+1)
    end if
  else if (k1 <= N1/2+1) then
    re =  x(1+2*(k1-1)+0, k2)
    im =  x(1+2*(k1-1)+1, k2)
  else
    re =  x(1+2*(N1-k1+1)+0, 1+N2-k2+1)
    im = -x(1+2*(N1-k1+1)+1, 1+N2-k2+1)
  end if

  res = cmplx(re, im, WP)
end subroutine ccs

! Fetch x(k1,k2) from the result of N1-by-N2 real FFT PACK-packed in x(:,:)
! Assume k1=1:N1, k2=1:N2
subroutine pack(res, x, LD1, N1,N2, k1,k2)

  complex(WP) :: res
  integer :: LD1, N1, N2, k1, k2
  real(WP) :: x(LD1,*)

  real(WP) :: re, im

  if (k2 == 1) then
    if (k1 == 1) then
      re =  x(1,1)
      im =  0
    else if (k1-1 == N1-k1+1) then
      re =  x(2*(k1-1),1)
      im =  0
    else if (k1 <= N1/2+1) then
      re =  x(2*(k1-1)+0,1)
      im =  x(2*(k1-1)+1,1)
    else
      re =  x(2*(N1-k1+1)+0,1)
      im = -x(2*(N1-k1+1)+1,1)
    end if
  else if (k1 == 1) then
    if (k2-1 == N2-k2+1) then
      re =  x(1,N2) !?
      im =  0
    else if (k2 <= N2/2+1) then
      re =  x(1,2*(k2-1)+0)
      im =  x(1,2*(k2-1)+1)
    else
      re =  x(1,2*(N2-k2+1)+0)
      im = -x(1,2*(N2-k2+1)+1)
    endif
  else if (k1-1 == N1-k1+1) then
    if (k2-1 == N2-k2+1) then
      re =  x(N1,N2)
      im =  0
    else if (k2 <= N2/2+1) then
      re =  x(N1,2*(k2-1)+0)
      im =  x(N1,2*(k2-1)+1)
    else
      re =  x(N1,2*(N2-k2+1)+0)
      im = -x(N1,2*(N2-k2+1)+1)
    end if
  else if (k1 <= N1/2+1) then
    re =  x(2*(k1-1)+0,1+k2-1)
    im =  x(2*(k1-1)+1,1+k2-1)
  else
    re =  x(2*(N1-k1+1)+0,1+N2-k2+1)
    im = -x(2*(N1-k1+1)+1,1+N2-k2+1)
  end if

  res = cmplx(re, im, WP)
end subroutine pack

! Fetch x(k1,k2) from the result of N1-by-N2 real FFT PERM-packed in x(:,:)
! Assume k1=1:N1, k2=1:N2
subroutine perm(res, x, LD1, N1,N2, k1,k2)

  complex(WP) :: res
  integer :: LD1, N1, N2, k1, k2
  real(WP) :: x(LD1, *)

  real(WP) :: re, im

  if (k2 == 1) then
    if (k1 == 1) then
      re =  x(1,1)
      im =  0
    else if (k1-1 == N1-k1+1) then
      re =  x(2,1)
      im =  0
    else if (k1 <= N1/2+1) then
      re =  x(1+2*(k1-1)+0 - mod(N1,2),1)
      im =  x(1+2*(k1-1)+1 - mod(N1,2),1)
    else
      re =  x(1+2*(N1-k1+1)+0 - mod(N1,2),1)
      im = -x(1+2*(N1-k1+1)+1 - mod(N1,2),1)
    end if
  else if (k1 == 1) then
    if (k2-1 == N2-k2+1) then
      re =  x(1,2)
      im =  0
    else if (k2 <= N2/2+1) then
      re =  x(1,1+2*(k2-1)+0 - mod(N2,2))
      im =  x(1,1+2*(k2-1)+1 - mod(N2,2))
    else
      re =  x(1,1+2*(N2-k2+1)+0 - mod(N2,2))
      im = -x(1,1+2*(N2-k2+1)+1 - mod(N2,2))
    endif
  else if (k1-1 == N1-k1+1) then
    if (k2-1 == N2-k2+1) then
      re =  x(2,2)
      im =  0
    else if (k2 <= N2/2+1) then
      re =  x(2,1+2*(k2-1)+0-mod(N2,2))
      im =  x(2,1+2*(k2-1)+1-mod(N2,2))
    else
      re =  x(2,1+2*(N2-k2+1)+0-mod(N2,2))
      im = -x(2,1+2*(N2-k2+1)+1-mod(N2,2))
    end if
  else if (k1 <= N1/2+1) then
    re =  x(1+2*(k1-1)+0-mod(N1,2),1+k2-1)
    im =  x(1+2*(k1-1)+1-mod(N1,2),1+k2-1)
  else
    re =  x(1+2*(N1-k1+1)+0-mod(N1,2),1+N2-k2+1)
    im = -x(1+2*(N1-k1+1)+1-mod(N1,2),1+N2-k2+1)
  end if

  res = cmplx(re, im, WP)
end subroutine perm


end module config_conjugate_even_storage
