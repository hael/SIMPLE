!@descr: a 2D complex field with its own c2c FFTW plans and the native-disk embedding, for wave-optics image formation
! (doc/implementation_notes/planned/ghostbuster_operators.md)
module simple_wave_field
use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_float_complex, c_ptr, c_size_t, c_null_ptr
use simple_core_module_api, only: simple_exception
use simple_fftw3,           only: fftwf_alloc_complex, fftwf_free, fftwf_plan_dft_2d, fftwf_execute_dft, &
    &fftwf_destroy_plan, FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
implicit none
private
#include "simple_local_flags.inc"

public :: wave_field, wave_freq

!> A complex field on an n x n grid, x the first index, in real space or Fourier space as the caller
!! keeps track; Fourier arrays are in FFT order (index i holds frequency wave_freq(i,n)). The transforms
!! are in place on the field's own FFTW buffer, with plans made once (FFTW_ESTIMATE, so the arithmetic is
!! the same in every run): fwd is the unnormalised forward DFT F (sign -1), bwd its inverse F^-1 = F^H/n^2,
!! so the transposes are fwd^H = n^2 bwd and bwd^H = fwd/n^2. One field belongs to one thread; its plans
!! are made and destroyed in the unnamed critical section, the lock every FFTW planner call of SIMPLE
!! (image) takes. A field is never copied: a copy would alias the FFTW buffer.
!!
!! A native-disk spectrum S(h,k), -R <= h,k <= R, h^2+k^2 <= R^2 (both signs of h, the plane of
!! pcg_lattice) is the transform of a native image of box nat = 2R. embed places it on the fine grid,
!! scaled by (n/nat)^2 so that bwd gives the native image's band-limited interpolant at its own
!! amplitude; crop takes the disk back, scaled by (nat/n)^2, so that crop(embed(S)) = S and crop of
!! fwd(f) is the native transform of f sampled at the native pixels. embed_t and crop_t are their
!! transposes under Re<a,b>. On a grid of the native box itself (n = nat) the disk's two Nyquist points of
!! an axis fall on one bin: embed adds them (the aliasing of the coarse grid) and crop gives the bin to
!! both, so the pairs stay transposes; crop(embed(S)) = S then holds off the Nyquist points.
type :: wave_field
    private
    integer                           :: n        = 0
    type(c_ptr)                       :: data_ptr = c_null_ptr
    type(c_ptr)                       :: plan_fwd = c_null_ptr
    type(c_ptr)                       :: plan_bwd = c_null_ptr
    complex(c_float_complex), pointer :: cmat(:,:) => null()
    complex(c_float_complex), pointer :: cout(:,:) => null() !< the same buffer: FFTW's output of an in-place plan
contains
    procedure :: new
    procedure :: get_n
    procedure :: ptr
    procedure :: zero
    procedure :: fwd
    procedure :: bwd
    procedure :: embed
    procedure :: embed_t
    procedure :: crop
    procedure :: crop_t
    procedure, private :: put_disk
    procedure, private :: get_disk
    procedure :: kill
end type wave_field

contains

    !> the signed frequency of FFT-order index i (1-based) of an n-point transform
    elemental integer function wave_freq( i, n )
        integer, intent(in) :: i, n
        wave_freq = i - 1
        if( wave_freq > n / 2 ) wave_freq = wave_freq - n
    end function wave_freq

    subroutine new( self, n )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: n
        if( n < 2 ) THROW_HARD('a wave field needs at least two points per side; new')
        call self%kill
        self%n        = n
        self%data_ptr = fftwf_alloc_complex(int(n, c_size_t) * int(n, c_size_t))
        if( .not. c_associated(self%data_ptr) ) THROW_HARD('FFTW allocation failed; wave_field%new')
        call c_f_pointer(self%data_ptr, self%cmat, [n, n])
        call c_f_pointer(self%data_ptr, self%cout, [n, n])
        self%cmat = cmplx(0., 0., kind=c_float_complex)
        !$omp critical
        self%plan_fwd = fftwf_plan_dft_2d(n, n, self%cmat, self%cout, FFTW_FORWARD,  FFTW_ESTIMATE)
        self%plan_bwd = fftwf_plan_dft_2d(n, n, self%cmat, self%cout, FFTW_BACKWARD, FFTW_ESTIMATE)
        !$omp end critical
        if( .not. c_associated(self%plan_fwd) .or. .not. c_associated(self%plan_bwd) )then
            THROW_HARD('FFTW plan creation failed; wave_field%new')
        endif
    end subroutine new

    pure integer function get_n( self )
        class(wave_field), intent(in) :: self
        get_n = self%n
    end function get_n

    !> the field's values, for the pointwise operations of its owner
    function ptr( self ) result( p )
        class(wave_field), intent(in) :: self
        complex(c_float_complex), pointer :: p(:,:)
        p => self%cmat
    end function ptr

    subroutine zero( self )
        class(wave_field), intent(inout) :: self
        self%cmat = cmplx(0., 0., kind=c_float_complex)
    end subroutine zero

    !> in place: the forward DFT, unnormalised
    subroutine fwd( self )
        class(wave_field), intent(inout) :: self
        call fftwf_execute_dft(self%plan_fwd, self%cmat, self%cout)
    end subroutine fwd

    !> in place: the inverse DFT, normalised by 1/n^2
    subroutine bwd( self )
        class(wave_field), intent(inout) :: self
        call fftwf_execute_dft(self%plan_bwd, self%cmat, self%cout)
        self%cmat = self%cmat * (1. / real(self%n)**2)
    end subroutine bwd

    !> the field (Fourier space) becomes the native-disk spectrum S, scaled by (n/nat)^2, zero elsewhere
    subroutine embed( self, s, R )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: R
        complex,           intent(in)    :: s(-R:R,-R:R)
        call self%put_disk(s, R, (real(self%n) / real(2*R))**2)
    end subroutine embed

    !> the transpose of embed: the field's disk, scaled by (n/nat)^2
    subroutine embed_t( self, s, R )
        class(wave_field), intent(in)  :: self
        integer,           intent(in)  :: R
        complex,           intent(out) :: s(-R:R,-R:R)
        call self%get_disk(s, R, (real(self%n) / real(2*R))**2)
    end subroutine embed_t

    !> the band-limited crop of the field (Fourier space) to the native disk, scaled by (nat/n)^2
    subroutine crop( self, s, R )
        class(wave_field), intent(in)  :: self
        integer,           intent(in)  :: R
        complex,           intent(out) :: s(-R:R,-R:R)
        call self%get_disk(s, R, (real(2*R) / real(self%n))**2)
    end subroutine crop

    !> the transpose of crop: the disk S on the field, scaled by (nat/n)^2, zero elsewhere
    subroutine crop_t( self, s, R )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: R
        complex,           intent(in)    :: s(-R:R,-R:R)
        call self%put_disk(s, R, (real(2*R) / real(self%n))**2)
    end subroutine crop_t

    subroutine put_disk( self, s, R, scale )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: R
        complex,           intent(in)    :: s(-R:R,-R:R)
        real,              intent(in)    :: scale
        integer :: h, k, i, j
        if( 2*R > self%n ) THROW_HARD('the native disk does not fit the field; put_disk')
        self%cmat = cmplx(0., 0., kind=c_float_complex)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) cycle
                i = modulo(h,self%n) + 1
                j = modulo(k,self%n) + 1
                self%cmat(i,j) = self%cmat(i,j) + scale * s(h,k)
            end do
        end do
    end subroutine put_disk

    subroutine get_disk( self, s, R, scale )
        class(wave_field), intent(in)  :: self
        integer,           intent(in)  :: R
        complex,           intent(out) :: s(-R:R,-R:R)
        real,              intent(in)  :: scale
        integer :: h, k
        if( 2*R > self%n ) THROW_HARD('the native disk does not fit the field; get_disk')
        s = cmplx(0., 0.)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) cycle
                s(h,k) = scale * self%cmat(modulo(h,self%n)+1, modulo(k,self%n)+1)
            end do
        end do
    end subroutine get_disk

    !> idempotent
    subroutine kill( self )
        class(wave_field), intent(inout) :: self
        if( c_associated(self%plan_fwd) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_fwd)
            !$omp end critical
        endif
        if( c_associated(self%plan_bwd) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_bwd)
            !$omp end critical
        endif
        if( c_associated(self%data_ptr) ) call fftwf_free(self%data_ptr)
        nullify(self%cmat, self%cout)
        self%data_ptr = c_null_ptr
        self%plan_fwd = c_null_ptr
        self%plan_bwd = c_null_ptr
        self%n        = 0
    end subroutine kill

end module simple_wave_field
