!@descr: utilities for convolution interpolation (gridding)
module simple_gridding
use simple_core_module_api
use simple_image,     only: image
use simple_projector, only: projector
implicit none

public :: prep2D_inv_instrfun4mul, prep3D_inv_kbenvelope4mul
public :: kb_stencil_envelope_1d, kb_stencil_inv_envelope_1d
public :: kb_stencil_centered_crop_inv_envelope_1d, deapodize3D_inplace
private
#include "simple_local_flags.inc"

contains

    !=============================================================
    ! 2D: SIMPLE convention => ldim_crop(3)=1, square (n x n)
    ! Returns inverse instrument function for multiplication.
    !=============================================================
    function prep2D_inv_instrfun4mul( ldim_crop, ldim_croppd, smpd_crop ) result( img )
        integer,      intent(in) :: ldim_crop(3), ldim_croppd(3)
        real,         intent(in) :: smpd_crop
        type(image)              :: img
        real(c_float), parameter :: EPS_DIV = 1.0e-8_c_float
        type(kbinterpol)         :: kbwin
        real(c_float)            :: center, dist, pad_sc, kbzero, kv
        integer                  :: i, j, n
        real(c_float)            :: inv_1d(ldim_crop(1))
        ! Determine single dimension (assumes square)
        n = ldim_crop(1)
        if (ldim_crop(2) /= n .or. ldim_crop(3) /= 1) then
            THROW_HARD('prep2D_inv_instrfun4mul: expected ldim_crop = [n,n,1]')
        end if
        ! Allocate target image to expected size
        call img%new(ldim_crop, smpd_crop)
        ! centre coordinate (Fortran 1-based)
        center = real(n/2 + 1, c_float)
        ! pad scaling: mirror original intent (use provided ldim_croppd(1))
        pad_sc = 1.0_c_float / real(ldim_croppd(1), c_float)
        ! create KB window object
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        ! normalization constant at zero (guard small values)
        kbzero = kbwin%instr(0.0_c_float)
        if ( abs(kbzero) < EPS_DIV ) kbzero = 1.0_c_float
        ! Precompute 1-D inverse vector
        do i = 1, n
            dist = pad_sc * ( real(i, c_float) - center )
            kv   = kbwin%instr(dist)
            if ( abs(kv) < EPS_DIV ) then
                inv_1d(i) = 0.0_c_float
            else
                inv_1d(i) = kbzero / kv
            end if
        end do
        ! Fill image with outer product inv_1d(i) * inv_1d(j)
        do j = 1, n
            do i = 1, n
                call img%set([i, j, 1], inv_1d(i) * inv_1d(j))
            end do
        end do
    end function prep2D_inv_instrfun4mul

    !=============================================================
    ! 3D gridding deapodization
    !
    ! The 3D reconstructor inserts on the NATIVE lattice (period box_crop)
    ! with the discrete, per-axis normalized KB stencil (half-width KBWINSZ
    ! native voxels). The real-space envelope multiplying the reconstructed
    ! map is therefore the inverse transform of that normalized stencil with
    ! period box_crop — not the continuous instrument function, and not a
    ! padded-lattice period (doc/implementation_notes/drop_legacy_box_division.md
    ! S2.1). The routines below are shared with the PCG backend, which uses
    ! the same stencil on its padded lattice (period box_croppd).
    !=============================================================

    !> 1-D real-space envelope of the discrete, per-axis normalized KB stencil
    !! at the Fourier origin (fractional offset 0) for a lattice of period n,
    !! normalized to unity at the centre (n/2+1).
    subroutine kb_stencil_envelope_1d( kbwin, n, env1d )
        type(kbinterpol),  intent(in)  :: kbwin
        integer,           intent(in)  :: n
        real, allocatable, intent(out) :: env1d(:)
        real, allocatable :: w(:)
        real    :: s, arg, ctrval
        integer :: wdim, iwinsz, i, q, x, c
        wdim   = kbwin%get_wdim()
        iwinsz = ceiling(kbwin%get_winsz() - 0.5)
        allocate(w(wdim))
        do i = 1, wdim
            w(i) = kbwin%apod(real(i - iwinsz - 1))
        end do
        s = sum(w)
        if( abs(s) > TINY ) w = w / s
        if( allocated(env1d) ) deallocate(env1d)
        allocate(env1d(n), source=0.)
        c = n/2 + 1
        do x = 1, n
            do i = 1, wdim
                q   = i - iwinsz - 1
                arg = 2.0 * PI * real(q * (x-c)) / real(n)
                env1d(x) = env1d(x) + w(i) * cos(arg)
            end do
        end do
        ctrval = env1d(c)
        if( abs(ctrval) > TINY ) env1d = env1d / ctrval
        deallocate(w)
    end subroutine kb_stencil_envelope_1d

    !> Guarded reciprocal of kb_stencil_envelope_1d for the production KB window
    subroutine kb_stencil_inv_envelope_1d( n, inv1d )
        integer,           intent(in)  :: n
        real, allocatable, intent(out) :: inv1d(:)
        real, parameter   :: EPS_DIV = 1.0e-8
        type(kbinterpol)  :: kbwin
        real, allocatable :: env1d(:)
        integer :: i
        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        call kb_stencil_envelope_1d(kbwin, n, env1d)
        if( allocated(inv1d) ) deallocate(inv1d)
        allocate(inv1d(n), source=0.)
        do i = 1, n
            if( abs(env1d(i)) >= EPS_DIV ) inv1d(i) = 1.0 / env1d(i)
        end do
        deallocate(env1d)
    end subroutine kb_stencil_inv_envelope_1d

    !> Guarded reciprocal of the padded-period KB envelope at the centered
    !! native-volume indices. The returned factors have native_box elements.
    subroutine kb_stencil_centered_crop_inv_envelope_1d( kbwin, padded_box, native_box, inv1d )
        type(kbinterpol),  intent(in)  :: kbwin
        integer,           intent(in)  :: padded_box, native_box
        real, allocatable, intent(out) :: inv1d(:)
        real, parameter :: EPS_DIV = 1.0e-8
        real, allocatable :: padded_env(:)
        integer :: i, offset
        if( native_box < 1 .or. native_box > padded_box .or. &
            &mod(padded_box-native_box,2) /= 0 ) &
            &error stop 'centered-crop envelope requires an even nonnegative size difference'
        call kb_stencil_envelope_1d(kbwin,padded_box,padded_env)
        offset = (padded_box-native_box)/2
        if( allocated(inv1d) ) deallocate(inv1d)
        allocate(inv1d(native_box),source=0.0)
        do i = 1, native_box
            if( abs(padded_env(offset+i)) >= EPS_DIV ) inv1d(i) = 1.0/padded_env(offset+i)
        end do
        deallocate(padded_env)
    end subroutine kb_stencil_centered_crop_inv_envelope_1d

    !> Multiplies a real-space cubic volume in place by the separable product
    !! inv1d(i)*inv1d(j)*inv1d(k)
    subroutine deapodize3D_inplace( vol, inv1d )
        class(image), intent(inout) :: vol
        real,         intent(in)    :: inv1d(:)
        real, pointer :: rmat(:,:,:) => null()
        integer :: ldim(3), i, j, k, n
        ldim = vol%get_ldim()
        n    = ldim(1)
        if( ldim(2) /= n .or. ldim(3) /= n ) THROW_HARD('deapodize3D_inplace: volume must be cubic')
        if( size(inv1d) /= n ) THROW_HARD('deapodize3D_inplace: envelope/volume size mismatch')
        if( vol%is_ft() ) THROW_HARD('deapodize3D_inplace: volume must be in real space')
        call vol%get_rmat_ptr(rmat)
        !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    rmat(i,j,k) = rmat(i,j,k) * (inv1d(i) * inv1d(j) * inv1d(k))
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rmat)
    end subroutine deapodize3D_inplace

    !> Inverse envelope volume (period ldim_crop(1)) for multiplication, for
    !! callers that finalize a bare reconstructor themselves (flex paths)
    function prep3D_inv_kbenvelope4mul( ldim_crop, smpd_crop ) result( img )
        integer, intent(in) :: ldim_crop(3)
        real,    intent(in) :: smpd_crop
        type(image)         :: img
        real, allocatable   :: inv1d(:)
        integer :: i, j, k, n
        n = ldim_crop(1)
        if( ldim_crop(2) /= n .or. ldim_crop(3) /= n )then
            THROW_HARD('prep3D_inv_kbenvelope4mul: ldim_crop must be cube (n x n x n)')
        end if
        call kb_stencil_inv_envelope_1d(n, inv1d)
        call img%new(ldim_crop, smpd_crop)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    call img%set([i, j, k], inv1d(i) * inv1d(j) * inv1d(k))
                end do
            end do
        end do
        deallocate(inv1d)
    end function prep3D_inv_kbenvelope4mul

end module simple_gridding
