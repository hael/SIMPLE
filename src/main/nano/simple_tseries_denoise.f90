!@descr: image-array adapter for FLEX kernel denoising of full time-series frames
module simple_tseries_denoise
use simple_defs,  only: sp, dp
use simple_error, only: simple_exception
use simple_image, only: image
implicit none
private
public :: flex_denoise_imgarr
#include "simple_local_flags.inc"

contains

    !> Denoise full frames without alignment or masking; retain their grid and individual means.
    subroutine flex_denoise_imgarr( imgs, den_imgs, ncomp, halfwin, l_eo )
        use simple_flex_traj_denoise, only: flex_traj_denoise
        type(image),              intent(inout) :: imgs(:)
        type(image), allocatable, intent(out)   :: den_imgs(:)
        integer,                  intent(in)    :: ncomp, halfwin
        logical,                  intent(in)    :: l_eo
        integer, parameter :: SEED = 20261026, MIN_FRAMES = 12
        type(image) :: work
        complex(sp), allocatable :: y(:,:), yden(:,:)
        real(sp),    allocatable :: wq(:), means(:)
        real(dp),    allocatable :: ess(:), choice(:,:)
        integer,     allocatable :: hidx(:), kidx(:), shell(:)
        integer :: n, ldim(3), i, j, h, k, nx, ny, ncoeff, kband, khere, kmate
        real :: smpd, shell_scale
        complex(sp) :: coeff
        logical :: boundary, self_conjugate
        n = size(imgs)
        if( n == 0 ) THROW_HARD('empty frame array; flex_denoise_imgarr')
        if( ncomp < 1 ) THROW_HARD('ncomp must be positive; flex_denoise_imgarr')
        if( halfwin < 0 ) THROW_HARD('halfwin must be nonnegative; flex_denoise_imgarr')
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( ldim(3) /= 1 ) THROW_HARD('2D frames required; flex_denoise_imgarr')
        if( any(mod(ldim(1:2),2) /= 0) ) THROW_HARD('even frame dimensions required; flex_denoise_imgarr')
        do i = 1, n
            if( any(imgs(i)%get_ldim() /= ldim) ) THROW_HARD('frame dimensions differ; flex_denoise_imgarr')
            if( abs(imgs(i)%get_smpd()-smpd) > 1.e-4 ) THROW_HARD('frame sampling differs; flex_denoise_imgarr')
            if( imgs(i)%is_ft() ) THROW_HARD('real-space frames required; flex_denoise_imgarr')
        enddo
        allocate(den_imgs(n))
        do i = 1, n
            call den_imgs(i)%new(ldim, smpd, wthreads=.false.)
        enddo
        if( n < MIN_FRAMES )then
            do i = 1, n
                call den_imgs(i)%copy_fast(imgs(i))
            enddo
            return
        endif
        nx = ldim(1)
        ny = ldim(2)
        ncoeff = nx*ny/2 + 1
        allocate(hidx(ncoeff), kidx(ncoeff), shell(ncoeff), wq(ncoeff), means(n))
        shell_scale = real(min(nx,ny))
        j = 0
        do h = 0, nx/2
            boundary = h == 0 .or. h == nx/2
            do k = -ny/2, ny/2-1
                if( h == 0 .and. k == 0 ) cycle
                if( boundary .and. k < 0 .and. k /= -ny/2 ) cycle
                j = j + 1
                hidx(j) = h
                kidx(j) = k
                ! Use physical radial frequency on rectangular grids; corners remain in the output.
                shell(j) = max(1,nint(shell_scale*sqrt((real(h)/real(nx))**2+(real(k)/real(ny))**2)))
                self_conjugate = boundary .and. (k == 0 .or. k == -ny/2)
                wq(j) = merge(1.,2.,self_conjugate)
            enddo
        enddo
        if( j /= ncoeff ) THROW_HARD('invalid Fourier lattice; flex_denoise_imgarr')
        allocate(y(ncoeff,n), yden(ncoeff,n), ess(n), choice(3,n))
        call work%new(ldim, smpd, wthreads=.false.)
        do i = 1, n
            call work%copy_fast(imgs(i))
            means(i) = work%mean()
            call work%add(-means(i))
            call work%fft
            do j = 1, ncoeff
                khere = modulo(kidx(j),ny) + 1
                y(j,i) = work%get_cmat_at(hidx(j)+1,khere,1)
            enddo
        enddo
        call work%kill
        ! The unmasked full frame supplies both the similarity and the kernel pre-image.
        call flex_traj_denoise(y, y, hidx, kidx, shell, wq, ncomp, 0, halfwin, l_eo, SEED, 1, &
            &yden, ess, choice, kband)
        do i = 1, n
            call den_imgs(i)%zero_and_flag_ft
            do j = 1, ncoeff
                h = hidx(j)
                k = kidx(j)
                khere = modulo(k,ny) + 1
                boundary = h == 0 .or. h == nx/2
                self_conjugate = boundary .and. (k == 0 .or. k == -ny/2)
                coeff = yden(j,i)
                if( self_conjugate ) coeff = cmplx(real(coeff,sp),0.,kind=sp)
                call den_imgs(i)%set_cmat_at(h+1,khere,1,coeff)
                if( boundary .and. .not. self_conjugate )then
                    kmate = modulo(-k,ny) + 1
                    call den_imgs(i)%set_cmat_at(h+1,kmate,1,conjg(coeff))
                endif
            enddo
            call den_imgs(i)%ifft
            call den_imgs(i)%add(means(i))
        enddo
    end subroutine flex_denoise_imgarr

end module simple_tseries_denoise
