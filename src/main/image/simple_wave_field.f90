!@descr: a 2D complex field with its own FFTW plans (c2c, and r2c/c2r for an optional real field) and the native-disk embedding
! (conventions: doc/policies/3D/wave_model_policy.md, section 3.1)
module simple_wave_field
use, intrinsic :: iso_c_binding, only: c_associated, c_f_pointer, c_float, c_float_complex, c_ptr, c_size_t, c_null_ptr
use simple_error,           only: simple_exception
use simple_fftw3,           only: fftwf_alloc_complex, fftwf_alloc_real, fftwf_free, fftwf_plan_dft_2d, fftwf_execute_dft, &
    &fftwf_plan_dft_r2c_2d, fftwf_plan_dft_c2r_2d, fftwf_execute_dft_r2c, fftwf_execute_dft_c2r, &
    &fftwf_destroy_plan, FFTW_FORWARD, FFTW_BACKWARD, FFTW_ESTIMATE
implicit none
private
#include "simple_local_flags.inc"

public :: wave_field, wave_freq

!> An n x n complex field (x the first index), in real or Fourier space as the caller keeps track, with in-place
!! transforms on its own FFTW buffer and the native-disk embed/crop and their transposes; one field per thread, never
!! copied (a copy would alias the buffer). Normalisations and scales: wave_model_policy.md, section 3.1
type :: wave_field
    private
    integer                           :: n        = 0
    type(c_ptr)                       :: data_ptr = c_null_ptr
    type(c_ptr)                       :: plan_fwd = c_null_ptr
    type(c_ptr)                       :: plan_bwd = c_null_ptr
    complex(c_float_complex), pointer :: cmat(:,:) => null()
    complex(c_float_complex), pointer :: cout(:,:) => null() !< the same buffer: FFTW's output of an in-place plan
    ! an optional real field and its half spectrum, out of place (new with_real)
    type(c_ptr)                       :: rdata_ptr = c_null_ptr, hdata_ptr = c_null_ptr
    type(c_ptr)                       :: plan_r2c  = c_null_ptr, plan_c2r  = c_null_ptr
    real(c_float),            pointer :: rmat(:,:) => null()  !< (n,n)
    complex(c_float_complex), pointer :: hmat(:,:) => null()  !< (0:n/2,n), the half spectrum in FFT order along k
contains
    procedure :: new
    procedure :: get_n
    procedure :: ptr
    procedure :: rptr
    procedure :: fwd_real_crop
    procedure :: crop_t_real
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

    !> with_real adds the real field of rptr, fwd_real_crop and crop_t_real
    subroutine new( self, n, with_real )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: n
        logical, optional, intent(in)    :: with_real
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
        if( .not. present(with_real) ) return
        if( .not. with_real ) return
        if( mod(n,2) /= 0 ) THROW_HARD('the real field needs an even n; wave_field%new')
        self%rdata_ptr = fftwf_alloc_real(int(n, c_size_t) * int(n, c_size_t))
        self%hdata_ptr = fftwf_alloc_complex(int(n/2+1, c_size_t) * int(n, c_size_t))
        if( .not. c_associated(self%rdata_ptr) .or. .not. c_associated(self%hdata_ptr) )then
            THROW_HARD('FFTW allocation failed; wave_field%new')
        endif
        call c_f_pointer(self%rdata_ptr, self%rmat, [n, n])
        call c_f_pointer(self%hdata_ptr, self%hmat, [n/2+1, n])
        self%hmat(0:,1:) => self%hmat
        self%rmat = 0.
        self%hmat = cmplx(0., 0., kind=c_float_complex)
        !$omp critical
        self%plan_r2c = fftwf_plan_dft_r2c_2d(n, n, self%rmat, self%hmat, FFTW_ESTIMATE)
        self%plan_c2r = fftwf_plan_dft_c2r_2d(n, n, self%hmat, self%rmat, FFTW_ESTIMATE)
        !$omp end critical
        if( .not. c_associated(self%plan_r2c) .or. .not. c_associated(self%plan_c2r) )then
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

    !> the real field, for the pointwise operations of its owner
    function rptr( self ) result( p )
        class(wave_field), intent(in) :: self
        real(c_float), pointer :: p(:,:)
        p => self%rmat
    end function rptr

    !> crop(fwd(f)) of the real field f: its forward DFT by r2c, the negative h from Hermitian symmetry
    subroutine fwd_real_crop( self, s, R )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: R
        complex,           intent(out)   :: s(-R:R,-R:R)
        real    :: scale
        integer :: h, k
        if( .not. c_associated(self%plan_r2c) ) THROW_HARD('no real field; wave_field%fwd_real_crop')
        if( 2*R > self%n ) THROW_HARD('the native disk does not fit the field; fwd_real_crop')
        call fftwf_execute_dft_r2c(self%plan_r2c, self%rmat, self%hmat)
        scale = (real(2*R) / real(self%n))**2
        s = cmplx(0., 0.)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) cycle
                if( h >= 0 )then
                    s(h,k) = scale * self%hmat(h, modulo(k,self%n)+1)
                else
                    s(h,k) = scale * conjg(self%hmat(-h, modulo(-k,self%n)+1))
                endif
            end do
        end do
    end subroutine fwd_real_crop

    !> the real field becomes n^2 Re(bwd(crop_t(s))): the Hermitian part of crop_t(s), inverted by c2r. Each
    !! disk point adds half its value to its bin and half its conjugate to the mirrored bin, wherever these
    !! fall on the half grid, so the Nyquist points that share a bin (n = nat) add as in crop_t
    subroutine crop_t_real( self, s, R )
        class(wave_field), intent(inout) :: self
        integer,           intent(in)    :: R
        complex,           intent(in)    :: s(-R:R,-R:R)
        real    :: scale
        integer :: h, k, hb, kb
        if( .not. c_associated(self%plan_c2r) ) THROW_HARD('no real field; wave_field%crop_t_real')
        if( 2*R > self%n ) THROW_HARD('the native disk does not fit the field; crop_t_real')
        scale = 0.5 * (real(2*R) / real(self%n))**2
        self%hmat = cmplx(0., 0., kind=c_float_complex)
        do k = -R, R
            do h = -R, R
                if( h*h + k*k > R*R ) cycle
                hb = modulo(h, self%n)
                kb = modulo(k, self%n)
                if( hb <= self%n/2 ) self%hmat(hb, kb+1) = self%hmat(hb, kb+1) + scale * s(h,k)
                hb = modulo(-h, self%n)
                kb = modulo(-k, self%n)
                if( hb <= self%n/2 ) self%hmat(hb, kb+1) = self%hmat(hb, kb+1) + scale * conjg(s(h,k))
            end do
        end do
        call fftwf_execute_dft_c2r(self%plan_c2r, self%hmat, self%rmat)
    end subroutine crop_t_real

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
        if( c_associated(self%plan_r2c) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_r2c)
            !$omp end critical
        endif
        if( c_associated(self%plan_c2r) )then
            !$omp critical
            call fftwf_destroy_plan(self%plan_c2r)
            !$omp end critical
        endif
        if( c_associated(self%rdata_ptr) ) call fftwf_free(self%rdata_ptr)
        if( c_associated(self%hdata_ptr) ) call fftwf_free(self%hdata_ptr)
        nullify(self%rmat, self%hmat)
        self%rdata_ptr = c_null_ptr
        self%hdata_ptr = c_null_ptr
        self%plan_r2c  = c_null_ptr
        self%plan_c2r  = c_null_ptr
        if( c_associated(self%data_ptr) ) call fftwf_free(self%data_ptr)
        nullify(self%cmat, self%cout)
        self%data_ptr = c_null_ptr
        self%plan_fwd = c_null_ptr
        self%plan_bwd = c_null_ptr
        self%n        = 0
    end subroutine kill

end module simple_wave_field
