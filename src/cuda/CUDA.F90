module CUDA
use cuda_h              ! include all
use cuda_runtime_h      ! include all
implicit none
contains

  subroutine my_cudaErrorCheck( err, b )
      use, intrinsic :: ISO_C_BINDING
      integer (KIND(cudaSuccess)), intent(in) :: err

    character, pointer :: string_err(:)
    integer :: i
    logical :: b

    print *, 'checking ...', err
    b = .false.
    if ( err /= cudaSuccess ) then
       print *, 'Error found'
       call C_F_POINTER( cudaGetErrorString( err ), string_err, [1] )
       i = 1
       do while ( string_err( i ) /= C_NULL_CHAR )
          i = i + 1
       enddo

       print *, 'CUDA Error Detected'
       print *, string_err( 1:i )
       b = .true.
    endif

  end subroutine my_cudaErrorCheck

end module
