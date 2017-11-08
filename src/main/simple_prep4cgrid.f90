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
    integer           :: lims(3,2)
    integer           :: recl
    type(image)       :: mskimg, img4grid
    logical           :: exists = .false.
contains
    procedure          :: new
    procedure          :: get_lims
    procedure          :: get_recl
    procedure, private :: memoize_instr_fun
    procedure, private :: prep_1
    procedure, private :: prep_2
    generic            :: prep => prep_1, prep_2
    procedure          :: kill
end type prep4cgrid

contains

    subroutine new( self, img, kbwin, ldim_pd )
        use simple_ftiter, only: ftiter
        class(prep4cgrid), intent(inout) :: self
        class(image),      intent(in)    :: img
        class(kbinterpol), intent(in)    :: kbwin
        integer,           intent(in)    :: ldim_pd(3)
        complex(sp), allocatable :: cmat(:,:)
        type(ftiter) :: fit 
        integer      :: npix
        call self%kill
        ! set logical dimensions of images
        self%ldim    = img%get_ldim()
        self%ldim_pd = ldim_pd
        ! checks
        if( .not. img%square_dims() ) stop 'ERROR, square dims assumed; simple_prep4cgrid :: new'
        if( .not. img%is_2d()       ) stop 'ERROR, only for 2D images; simple_prep4cgrid :: new'
        allocate( self%instr_fun(self%ldim_pd(1),self%ldim_pd(1),1) )
        ! memoize mask image and instrument function
        call self%mskimg%disc(self%ldim, img%get_smpd(), real(self%ldim(1)/2-2))
        call self%memoize_instr_fun(kbwin)
        ! work out size of complex matrix output
        call fit%new(ldim_pd, img%get_smpd())
        self%lims = fit%loop_lims(3)
        ! allocate complex matrix and check record lenght for writing to binary direct access file
        allocate(cmat(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2)))
        inquire(iolength=self%recl) cmat
        deallocate(cmat)
        ! make an img4grid for the case where the return value is a complex matrix
        call self%img4grid%new(self%ldim_pd, img%get_smpd())
        ! flag existence
        self%exists = .true.
    end subroutine new

    function get_lims( self ) result( lims )
        class(prep4cgrid), intent(in) :: self
        integer :: lims(3,2)
        lims = self%lims
    end function get_lims

    function get_recl( self ) result( recl )
        class(prep4cgrid), intent(in) :: self
        integer :: recl
        recl = self%recl
    end function get_recl

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
    subroutine prep_1( self, img, img4grid )
        class(prep4cgrid), intent(in)    :: self
        class(image),      intent(inout) :: img, img4grid
        call img%subtr_backgr_pad_divwinstr_fft(self%mskimg, self%instr_fun, img4grid)
    end subroutine prep_1

    !>  \brief  prepare cmat for gridding interpolation in Fourier space
    subroutine prep_2( self, img, cmat )
        class(prep4cgrid), intent(inout) :: self
        class(image),      intent(inout) :: img
        complex(sp),       intent(out)   :: cmat(self%lims(1,1):self%lims(1,2),self%lims(2,1):self%lims(2,2))
        integer :: h, k, logi(3), phys(3)
        call img%subtr_backgr_pad_divwinstr_fft(self%mskimg, self%instr_fun, self%img4grid)
        do h=self%lims(1,1),self%lims(1,2)
            do k=self%lims(2,1),self%lims(2,2)
                logi      = [h,k,0]
                phys      = self%img4grid%comp_addr_phys(logi)
                cmat(h,k) = self%img4grid%get_fcomp(logi,phys)
            end do
        end do
    end subroutine prep_2

    subroutine kill( self )
        class(prep4cgrid), intent(inout) :: self
        if( self%exists )then
            call self%mskimg%kill
            call self%prep4cgrid%kill
            deallocate(self%instr_fun)
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prep4cgrid