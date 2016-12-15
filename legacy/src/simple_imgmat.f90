!>  \brief  SIMPLE image matrix class
module simple_imgmat
use simple_fftw3 ! singleton
implicit none

public :: imgmat
private

type :: imgmat
    type(c_ptr), private                   :: p                   !< c pointer for fftw allocation
    real(kind=c_float), pointer            :: rmat(:,:,:)=>null() !< image pixels/voxels (in data)
    complex(kind=c_float_complex), pointer :: cmat(:,:,:)=>null() !< Fourier components
    logical                                :: exists=.false.      !< to indicate existence
  contains
    procedure :: allocate
    procedure :: deallocate
end type

contains

    subroutine allocate( self, ldim, which )
        use simple_math, only: fdim
        class(imgmat), intent(inout) :: self
        character(len=*), intent(in) :: which
        integer, intent(in)          :: ldim(3)
        integer                      :: array_shape(3)
        call self%deallocate
        array_shape(1)   = fdim(ldim(1))
        array_shape(2:3) = ldim(2:3)
        if( which .eq. 'complex' )then
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            self%p = fftwf_alloc_complex(int(product(array_shape),c_size_t))
            ! Set up the complex array which will point at the allocated memory
            call c_f_pointer(self%p,self%cmat,array_shape)
            self%cmat = cmplx(0.,0.)
        else if( which .eq. 'real' )then
            ! Letting FFTW do the allocation in C ensures that we will be using aligned memory
            self%p = fftwf_alloc_real(int(product(array_shape),c_size_t))
            ! Set up the real array
            call c_f_pointer(self%p,self%rmat,array_shape)
            self%rmat = 0.
        else
            stop 'unknown type (which) in simple_imgmat::allocate'
        endif
        self%exists = .true.
    end subroutine
    
    subroutine deallocate( self )
        class(imgmat), intent(inout) :: self
        if( self%exists )then
            call fftwf_free(self%p)
            self%rmat=>null()
            self%cmat=>null()
        endif
    end subroutine
    
end module
