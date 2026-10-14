!@descr: the memoised plane grid of gen_fplane4rec: the shell lookup table of a padded Fourier plane and, per
!  even/odd half, the noise spectrum resampled to its shells as an inverse and an inverse square root
! A reconstruction plane is built per particle on a padded grid fixed for the whole pass; building the
! shell table, resampling the spectrum and dividing by it per particle repeats the same work. The grid
! is built once (new), the spectra are attached once per spectra state (set_sigma2), and every plane of
! the pass reads both.
module simple_fplane_grid
use simple_defs,    only: OSMPL_PAD_FAC
use simple_error,   only: simple_exception
use simple_math_ft, only: resample_sigma2
implicit none
private
#include "simple_local_flags.inc"

public :: fplane_grid

type :: fplane_grid
    integer :: box_croppd = 0                     !< padded box of the plane grid
    integer :: nyq        = 0                     !< its Nyquist shell, box_croppd/2
    integer :: max_r2     = -1                    !< largest h*h+k*k the shell table covers
    integer :: kfrom      = 0                     !< lowest shell of the attached spectra (lower shells take its value)
    integer, allocatable :: shell_lut(:)          !< nint(sqrt(r2)) for r2 = 0:max_r2
    real,    allocatable :: inv_sig2(:,:)         !< 1/sigma2 at the padded shells 0:nyq, per half (column 1 even, 2 odd)
    real,    allocatable :: rsqrt_sig2(:,:)       !< 1/sqrt(sigma2), same layout
contains
    procedure :: new
    procedure :: set_sigma2
    procedure :: has_sigma2
    procedure :: matches
    procedure :: kill
end type fplane_grid

contains

    !> the grid of a padded plane of box box_croppd: its shell table, without spectra
    subroutine new( self, box_croppd )
        class(fplane_grid), intent(inout) :: self
        integer,            intent(in)    :: box_croppd
        integer :: r2
        call self%kill
        if( box_croppd < 2 ) THROW_HARD('plane grid box must be at least 2')
        self%box_croppd = box_croppd
        self%nyq        = box_croppd / 2
        ! covers every sample of the padded plane, |h|,|k| <= nyq+1
        self%max_r2     = 2 * (self%nyq + 1)**2
        allocate(self%shell_lut(0:self%max_r2))
        do r2 = 0, self%max_r2
            self%shell_lut(r2) = nint(sqrt(real(r2)))
        end do
    end subroutine new

    !> attaches the two spectra sigma2(kstart:,1:2) (native shells, column = half) as gen_fplane4rec resampled
    !! them: shells kfrom up to the native Nyquist box_croppd/(2*OSMPL_PAD_FAC), stretched over the padded
    !! shells 0:nyq, shells below kfrom taking the value at kfrom
    subroutine set_sigma2( self, kstart, sigma2, kfrom )
        class(fplane_grid), intent(inout) :: self
        integer,            intent(in)    :: kstart
        real,               intent(in)    :: sigma2(kstart:,:)
        integer,            intent(in)    :: kfrom
        real, allocatable :: resampled(:)
        character(len=:), allocatable :: msg
        integer :: sigma_nyq, ih
        if( .not. allocated(self%shell_lut) ) THROW_HARD('the plane grid has not been built; set_sigma2')
        if( size(sigma2,2) /= 2 ) THROW_HARD('one spectrum per even/odd half; set_sigma2')
        sigma_nyq = (self%box_croppd / OSMPL_PAD_FAC) / 2
        if( kfrom < kstart .or. sigma_nyq > ubound(sigma2,1) .or. kfrom > sigma_nyq )then
            msg = 'the spectra do not cover the shells of the plane grid; set_sigma2'
            THROW_HARD(msg)
        endif
        if( allocated(self%inv_sig2)   ) deallocate(self%inv_sig2)
        if( allocated(self%rsqrt_sig2) ) deallocate(self%rsqrt_sig2)
        allocate(self%inv_sig2(0:self%nyq,2), self%rsqrt_sig2(0:self%nyq,2), resampled(0:self%nyq))
        do ih = 1, 2
            call resample_sigma2(kfrom, sigma_nyq, sigma2(kfrom:sigma_nyq,ih), self%nyq, &
                &real(sigma_nyq)/real(self%nyq), resampled)
            self%inv_sig2(:,ih)   = 1.0 / resampled
            self%rsqrt_sig2(:,ih) = 1.0 / sqrt(resampled)
        end do
        self%kfrom = kfrom
        deallocate(resampled)
    end subroutine set_sigma2

    logical function has_sigma2( self )
        class(fplane_grid), intent(in) :: self
        has_sigma2 = allocated(self%inv_sig2)
    end function has_sigma2

    !> the grid of a padded box with the spectra attached from shell kfrom
    logical function matches( self, box_croppd, kfrom )
        class(fplane_grid), intent(in) :: self
        integer,            intent(in) :: box_croppd, kfrom
        matches = self%box_croppd == box_croppd .and. self%kfrom == kfrom .and. self%has_sigma2()
    end function matches

    subroutine kill( self )
        class(fplane_grid), intent(inout) :: self
        if( allocated(self%shell_lut)  ) deallocate(self%shell_lut)
        if( allocated(self%inv_sig2)   ) deallocate(self%inv_sig2)
        if( allocated(self%rsqrt_sig2) ) deallocate(self%rsqrt_sig2)
        self%box_croppd = 0
        self%nyq        = 0
        self%max_r2     = -1
        self%kfrom      = 0
    end subroutine kill

end module simple_fplane_grid
