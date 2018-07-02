module CUDA
use cuda_h              ! include all
use cuda_runtime_h      ! include all
implicit none
contains

    subroutine my_cudaErrorCheck( err, b , srcfile, srcline)
        use, intrinsic :: ISO_C_BINDING
        integer, intent(in) :: err
        character(len=*), intent(in), optional :: srcfile
        integer, intent(in), optional :: srcline
        character, pointer :: string_err(:)
        character(len=1024) :: sourcefile
        integer :: i, sourceline
        logical :: b
        integer (KIND(cudaSuccess)) :: errCuda
        errCuda = err
        b = .false.
        if ( err .ne. cudaSuccess ) then
            print *, 'Error found'
            call C_F_POINTER( cudaGetErrorString( err ), string_err, [1] )
            i = 1
            do while ( string_err( i ) /= C_NULL_CHAR )
                i = i + 1
            enddo

            print *, 'CUDA Error Detected'
            print *, string_err( 1:i )
            if(present(srcfile).and.present(srcline))then
                print *, 'CUDA error called from file:', trim(srcfile),':',srcline
            endif
            b = .true.
        endif


    end subroutine my_cudaErrorCheck

end module
