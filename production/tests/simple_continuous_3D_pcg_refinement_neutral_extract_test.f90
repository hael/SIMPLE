!@descr: Phase 2 neutral Cartesian Fourier and envelope extraction regression checks
module continuous_3D_pcg_refinement_neutral_extract_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true
use simple_cartesian_fourier, only: center_embed_real3d, center_crop_real3d, &
    &extract_native_fourier_plane, gather_packed_window, gather_packed_window_grad
use simple_defs, only: KBALPHA, KBWINSZ
use simple_gridding, only: kb_stencil_envelope_1d, &
    &kb_stencil_centered_crop_inv_envelope_1d
use simple_image, only: image
use simple_kbinterpol, only: kbinterpol
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_neutral_extract

integer, parameter :: BOX = 8
integer, parameter :: BOXPD = 16
integer, parameter :: WDIM = 3
real, parameter :: ENVELOPE_TOL = 2.e-6

contains

!> Compare each Phase 2 neutral operation with the pre-extraction algorithm.
subroutine run_neutral_extract()
    call test_centered_embed_crop()
    call test_centered_crop_envelope()
    call test_packed_gathers()
    call test_native_plane_extraction()
    write(*,'(a)') 'CONTINUOUS_3D_NEUTRAL_EXTRACT: PASS'
end subroutine run_neutral_extract

!> Verify exact centered indices and the embed/crop adjoint pair.
subroutine test_centered_embed_crop()
    real :: native(BOX,BOX,BOX), arbitrary_padded(BOXPD,BOXPD,BOXPD)
    real, allocatable :: legacy_crop(:,:,:), legacy_embed(:,:,:)
    real, allocatable :: neutral_crop(:,:,:), neutral_embed(:,:,:)
    integer :: i, j, k

    do k = 1, BOX
        do j = 1, BOX
            do i = 1, BOX
                native(i,j,k) = real(100*i+10*j+k)
            end do
        end do
    end do
    do k = 1, BOXPD
        do j = 1, BOXPD
            do i = 1, BOXPD
                arbitrary_padded(i,j,k) = real(10000*i+100*j+k)
            end do
        end do
    end do

    neutral_embed = center_embed_real3d(native,BOXPD)
    legacy_embed = legacy_center_embed(native,BOXPD)
    call assert_true(all(neutral_embed == legacy_embed), &
        &'neutral centered embed changed an index or value')
    neutral_crop = center_crop_real3d(neutral_embed,BOX)
    call assert_true(all(neutral_crop == native), &
        &'neutral crop did not invert the centered embed')
    neutral_crop = center_crop_real3d(arbitrary_padded,BOX)
    legacy_crop = legacy_center_crop(arbitrary_padded,BOX)
    call assert_true(all(neutral_crop == legacy_crop), &
        &'neutral centered crop changed an index or value')

    write(*,'(a)') 'CONTINUOUS_3D_NEUTRAL centered embed/crop: exact'
end subroutine test_centered_embed_crop

!> Verify padded-period native factors and the unchanged PCG 3-D envelope.
subroutine test_centered_crop_envelope()
    type(kbinterpol) :: kbwin
    type(reconstructor_pcg) :: pcgop
    real, allocatable :: env_padded(:), inv1d(:), legacy_inv1d(:)
    real, allocatable :: legacy_env3(:,:,:), legacy_invenv3(:,:,:)
    real, allocatable :: pcg_env(:,:,:), pcg_invenv(:,:,:)
    real :: center_value, max_env_error, max_factor_error, max_invenv_error
    integer :: i, j, k, offset

    kbwin = kbinterpol(KBWINSZ,KBALPHA)
    call kb_stencil_envelope_1d(kbwin,BOXPD,env_padded)
    call kb_stencil_centered_crop_inv_envelope_1d(kbwin,BOXPD,BOX,inv1d)
    offset = (BOXPD-BOX)/2
    allocate(legacy_inv1d(BOX),source=0.0)
    do i = 1, BOX
        if( abs(env_padded(offset+i)) >= 1.e-8 ) legacy_inv1d(i) = 1.0/env_padded(offset+i)
    end do
    max_factor_error = maxval(abs(inv1d-legacy_inv1d))
    call assert_true(max_factor_error == 0.0, &
        &'centered-crop inverse-envelope factors changed the old calculation')

    allocate(legacy_env3(BOX,BOX,BOX),legacy_invenv3(BOX,BOX,BOX))
    do k = 1, BOX
        do j = 1, BOX
            do i = 1, BOX
                legacy_env3(i,j,k) = env_padded(offset+i)*env_padded(offset+j)*env_padded(offset+k)
            end do
        end do
    end do
    center_value = legacy_env3(BOX/2+1,BOX/2+1,BOX/2+1)
    if( abs(center_value) < 1.e-8 ) center_value = 1.0
    legacy_env3 = legacy_env3/center_value
    legacy_invenv3 = 1.0
    where( abs(legacy_env3) < 1.e-8 )
        legacy_invenv3 = 0.0
    elsewhere
        legacy_invenv3 = 1.0/legacy_env3
    end where

    call pcgop%new(BOX,1.0)
    pcg_env = pcgop%get_env()
    pcg_invenv = pcgop%get_invenv()
    max_env_error = maxval(abs(pcg_env-legacy_env3))
    max_invenv_error = maxval(abs(pcg_invenv-legacy_invenv3)/max(1.0,abs(legacy_invenv3)))
    call assert_true(max_env_error <= ENVELOPE_TOL, &
        &'PCG envelope changed after the neutral extraction')
    call assert_true(max_invenv_error <= ENVELOPE_TOL, &
        &'PCG inverse envelope changed beyond single-precision roundoff')

    write(*,'(a,3(es14.6,1x))') 'CONTINUOUS_3D_NEUTRAL factor/env-abs/invenv-rel max: ', &
        &max_factor_error,max_env_error,max_invenv_error
    call pcgop%kill
end subroutine test_centered_crop_envelope

!> Verify the 27-tap order, packed/Friedel map, and negative wrap-table bound.
subroutine test_packed_gathers()
    complex :: cmat(BOXPD/2+1,BOXPD,BOXPD)
    complex :: legacy_gradient(3), legacy_value, neutral_gradient(3), neutral_value
    complex :: neutral_value_only
    real :: dw(WDIM,WDIM,WDIM,3), w(WDIM,WDIM,WDIM)
    integer :: i, j, k, axis, i0(3), wrap(-BOXPD-2:BOXPD+2)

    do k = 1, BOXPD
        do j = 1, BOXPD
            do i = 1, BOXPD/2+1
                cmat(i,j,k) = cmplx(real(100*i+3*j+k),real(-2*i+j-4*k))
            end do
        end do
    end do
    do i = lbound(wrap,1), ubound(wrap,1)
        wrap(i) = modulo(i+BOXPD/2,BOXPD)-BOXPD/2
    end do
    do k = 1, WDIM
        do j = 1, WDIM
            do i = 1, WDIM
                w(i,j,k) = real(i+2*j+3*k)/324.0
                do axis = 1, 3
                    dw(i,j,k,axis) = real((i-axis)*(j+axis)-k)/1000.0
                end do
            end do
        end do
    end do
    i0 = [-9,6,-2]

    neutral_value_only = gather_packed_window(cmat,lbound(wrap,1),wrap,i0,w)
    call gather_packed_window_grad(cmat,lbound(wrap,1),wrap,i0,w,dw, &
        &neutral_value,neutral_gradient)
    call legacy_packed_gather(cmat,lbound(wrap,1),wrap,i0,w,dw, &
        &legacy_value,legacy_gradient)
    call assert_true(neutral_value_only == legacy_value .and. neutral_value == legacy_value, &
        &'neutral packed value gather changed the old traversal')
    call assert_true(all(neutral_gradient == legacy_gradient), &
        &'neutral packed gradient gather changed the old traversal')

    write(*,'(a,i0,a,i0)') 'CONTINUOUS_3D_NEUTRAL packed gather wrap bounds: ', &
        &lbound(wrap,1),':',ubound(wrap,1)
end subroutine test_packed_gathers

!> Verify the neutral extraction and retained PCG wrapper against old code.
subroutine test_native_plane_extraction()
    type(image) :: img2d
    type(reconstructor_pcg) :: pcgop
    real :: pixels(BOX,BOX,1)
    complex, allocatable :: legacy_plane(:,:), neutral_plane(:,:), wrapper_plane(:,:)
    integer :: i, j, lims2(2,2), sqlp

    do j = 1, BOX
        do i = 1, BOX
            pixels(i,j,1) = real(3*i*i+5*j+2*i*j)
        end do
    end do
    call img2d%new([BOX,BOX,1],1.0)
    call img2d%set_rmat(pixels,.false.)
    call img2d%fft()
    call pcgop%new(BOX,1.0)
    lims2 = pcgop%get_lims2()
    sqlp = (BOX/2)**2
    neutral_plane = extract_native_fourier_plane(img2d,lims2,sqlp)
    wrapper_plane = pcgop%extract_native_plane(img2d)
    legacy_plane = legacy_native_plane(img2d,lims2,sqlp)
    call assert_true(all(neutral_plane == legacy_plane), &
        &'neutral native-plane extraction changed the old result')
    call assert_true(all(wrapper_plane == legacy_plane), &
        &'PCG native-plane wrapper changed the old result')

    write(*,'(a)') 'CONTINUOUS_3D_NEUTRAL native plane: exact'
    call pcgop%kill
    call img2d%kill
end subroutine test_native_plane_extraction

!> Pre-extraction centered embedding retained as a direct comparison oracle.
function legacy_center_embed( native, padded_box ) result( padded )
    real,    intent(in) :: native(:,:,:)
    integer, intent(in) :: padded_box
    real, allocatable :: padded(:,:,:)
    integer :: native_box, offset
    native_box = size(native,1)
    offset = (padded_box-native_box)/2
    allocate(padded(padded_box,padded_box,padded_box),source=0.0)
    padded(offset+1:offset+native_box,offset+1:offset+native_box, &
        &offset+1:offset+native_box) = native
end function legacy_center_embed

!> Pre-extraction centered crop retained as a direct comparison oracle.
function legacy_center_crop( padded, native_box ) result( native )
    real,    intent(in) :: padded(:,:,:)
    integer, intent(in) :: native_box
    real, allocatable :: native(:,:,:)
    integer :: offset
    offset = (size(padded,1)-native_box)/2
    allocate(native(native_box,native_box,native_box), &
        &source=padded(offset+1:offset+native_box,offset+1:offset+native_box, &
        &offset+1:offset+native_box))
end function legacy_center_crop

!> Pre-extraction packed gather retained as a direct comparison oracle.
pure subroutine legacy_packed_gather( cmat, wrap_lower, wrap, i0, w, dw, value, gradient )
    complex, intent(in)  :: cmat(:,:,:)
    integer, intent(in)  :: wrap_lower, wrap(wrap_lower:), i0(3)
    real,    intent(in)  :: w(:,:,:), dw(:,:,:,:)
    complex, intent(out) :: value, gradient(3)
    complex :: fcomp
    integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
    ny = size(cmat,2)
    nz = size(cmat,3)
    value = cmplx(0.,0.)
    gradient = cmplx(0.,0.)
    do dk = 1, size(w,3)
        mm = wrap(i0(3)+dk-1)
        do dj = 1, size(w,2)
            kk = wrap(i0(2)+dj-1)
            do di = 1, size(w,1)
                hh = wrap(i0(1)+di-1)
                if( hh >= 0 )then
                    ph = hh+1
                    pk = kk+1; if( kk < 0 ) pk = pk+ny
                    pm = mm+1; if( mm < 0 ) pm = pm+nz
                    fcomp = cmat(ph,pk,pm)
                else
                    ph = -hh+1
                    pk = -kk+1; if( -kk < 0 ) pk = pk+ny
                    pm = -mm+1; if( -mm < 0 ) pm = pm+nz
                    fcomp = conjg(cmat(ph,pk,pm))
                endif
                value = value+w(di,dj,dk)*fcomp
                gradient = gradient+dw(di,dj,dk,:)*fcomp
            end do
        end do
    end do
end subroutine legacy_packed_gather

!> Pre-extraction native-plane loop retained as a direct comparison oracle.
function legacy_native_plane( img2d, lims2, sqlp ) result( plane )
    class(image), intent(in) :: img2d
    integer,      intent(in) :: lims2(2,2), sqlp
    complex :: plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2))
    integer :: h, k, phys(3)
    plane = cmplx(0.,0.)
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k > sqlp ) cycle
            phys = img2d%comp_addr_phys(h,k,0)
            plane(h,k) = img2d%get_fcomp([h,k,0],phys)
        end do
    end do
end function legacy_native_plane

end module continuous_3D_pcg_refinement_neutral_extract_test
