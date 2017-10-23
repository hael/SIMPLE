module simple_prep4cgrid
#include "simple_lib.f08"
use simple_kbinterpol, only: kbinterpol
use simple_image,      only: image
implicit none

public :: prep4cgrid
private

type prep4cgrid
    private
    real, allocatable :: instr_fun(:,:,:)
    integer           :: ldim(3)
    integer           :: ldim_pd(3)
    type(image)       :: mskimg
    logical           :: exists = .false.
contains
    procedure          :: new
    procedure, private :: memoize_instr_fun
    procedure          :: prep
    procedure          :: kill
end type prep4cgrid

contains

    subroutine new( self, img, kbwin, ldim_pd )
        class(prep4cgrid), intent(inout) :: self
        class(image),      intent(in)    :: img
        class(kbinterpol), intent(in)    :: kbwin
        integer,           intent(in)    :: ldim_pd(3)
        integer :: npix
        call self%kill
        ! set logical dimensions of images
        self%ldim    = img%get_ldim()
        self%ldim_pd = ldim_pd
        ! checks
        if( .not. img%square_dims() ) stop 'ERROR, square dims assumed; simple_prep4cgrid :: new'
        if( .not. img%is_2d()       ) stop 'ERROR, only for 2D images; simple_prep4cgrid :: new'
        allocate( self%instr_fun(self%ldim_pd(1),self%ldim_pd(1),1) )
        call self%mskimg%disc(self%ldim, img%get_smpd(), real(self%ldim(1)/2-2))
        call self%memoize_instr_fun(kbwin)
        self%exists = .true.
    end subroutine new

    subroutine memoize_instr_fun( self, kbwin )
        class(prep4cgrid), intent(inout) :: self
        class(kbinterpol), intent(in)    :: kbwin
        real    :: w(self%ldim_pd(1))
        real    :: ci, arg
        integer :: i, j
        ci = -real(self%ldim_pd(1))/2.
        do i=1,self%ldim_pd(1)
            arg  = ci/real(self%ldim_pd(1))
            w(i) = kbwin%instr(arg)
            ci   = ci + 1.
        end do
        forall(i=1:self%ldim_pd(1), j=1:self%ldim_pd(1) ) self%instr_fun(i,j,1) = w(i)*w(j)
    end subroutine memoize_instr_fun

    !>  \brief  prepare image for gridding interpolation in Fourier space
    subroutine prep( self, img, img4grid )
        class(prep4cgrid), intent(in)    :: self
        class(image),      intent(inout) :: img, img4grid
        call img%subtr_backgr_pad_divwinstr_fft(self%mskimg, self%instr_fun, img4grid)
    end subroutine prep

    subroutine kill( self )
        class(prep4cgrid), intent(inout) :: self
        if( self%exists )then
            call self%mskimg%kill
            deallocate(self%instr_fun)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prep4cgrid