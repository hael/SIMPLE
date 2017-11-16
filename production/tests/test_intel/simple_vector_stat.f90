#ifndef INTEL

module simple_vector_stat
end module simple_vector_stat
#else
    include 'mkl_vsl.f90'

module simple_vector_stat
    use mkl_vsl
    implicit none
contains

    subroutine MKL_VSL_GAUSSIAN

      USE MKL_VSL_TYPE
      USE MKL_VSL

      real(kind=8) r(1000)  ! buffer for random numbers
      real(kind=8) s        ! average
      real(kind=8) a, sigma ! parameters of normal distribution

      TYPE (VSL_STREAM_STATE) :: stream

      integer(kind=4) errcode
      integer(kind=4) i,j
      integer brng,method,seed,n

      n = 1000
      s = 0.0
      a = 5.0
      sigma  = 2.0
      brng=VSL_BRNG_MT19937
      method=VSL_RNG_METHOD_GAUSSIAN_ICDF
      seed=777

!     ***** Initializing *****
      errcode=vslnewstream( stream, brng,  seed )

!     ***** Generating *****
      do i = 1,10
          errcode=vdrnggaussian( method, stream, n, r, a, sigma )
          do j = 1, 1000
              s = s + r(j)
          end do
      end do

      s = s / 10000.0

!     ***** Deinitialize *****
      errcode=vsldeletestream( stream )

!     ***** Printing results *****
      print *,"Sample mean of normal distribution = ", s

  end subroutine MKL_VSL_GAUSSIAN

end module simple_vector_stat
#endif
