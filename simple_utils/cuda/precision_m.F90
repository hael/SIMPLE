! 
!     Copyright (c) 2016, NVIDIA CORPORATION.  All rights reserved.
!
! NVIDIA CORPORATION and its licensors retain all intellectual property
! and proprietary rights in and to this software, related documentation
! and any modifications thereto.
!
!
!    These example codes are a portion of the code samples from the companion
!    website to the book "CUDA Fortran for Scientists and Engineers":
!
! http://store.elsevier.com/product.jsp?isbn=9780124169708
!

module precision_m
  integer, parameter :: singlePrecision = kind(0.0)
  integer, parameter :: doublePrecision = kind(0.0d0)
  
#ifdef DOUBLE
  integer, parameter :: fp_kind = doublePrecision
#else
  integer, parameter :: fp_kind = singlePrecision
#endif
end module precision_m
