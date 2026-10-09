!@descr: neutral Cartesian Fourier lattice embedding, extraction, and packed KB gathers
module simple_cartesian_fourier
use simple_image, only: image
use simple_kbinterpol, only: kbinterpol
implicit none
private

public :: center_embed_real3d, center_crop_real3d
public :: extract_native_fourier_plane
public :: gather_packed_window, gather_packed_window_grad, gather_packed_kb3_window_grad

contains

    !> Center a cubic real array in a larger cubic array. The size difference
    !! must be even so that the centered offset is an integer.
    function center_embed_real3d( native, padded_box ) result( padded )
        real,    intent(in) :: native(:,:,:)
        integer, intent(in) :: padded_box
        real, allocatable :: padded(:,:,:)
        integer :: native_box, offset
        native_box = size(native,1)
        if( size(native,2) /= native_box .or. size(native,3) /= native_box ) &
            &error stop 'center_embed_real3d requires a cubic input array'
        if( padded_box < native_box .or. mod(padded_box-native_box,2) /= 0 ) &
            &error stop 'center_embed_real3d requires an even nonnegative size difference'
        offset = (padded_box-native_box)/2
        allocate(padded(padded_box,padded_box,padded_box), source=0.0)
        padded(offset+1:offset+native_box, offset+1:offset+native_box, &
            &offset+1:offset+native_box) = native
    end function center_embed_real3d

    !> Return the centered native cube from a larger cubic real array.
    function center_crop_real3d( padded, native_box ) result( native )
        real,    intent(in) :: padded(:,:,:)
        integer, intent(in) :: native_box
        real, allocatable :: native(:,:,:)
        integer :: padded_box, offset
        padded_box = size(padded,1)
        if( size(padded,2) /= padded_box .or. size(padded,3) /= padded_box ) &
            &error stop 'center_crop_real3d requires a cubic input array'
        if( native_box < 1 .or. native_box > padded_box .or. mod(padded_box-native_box,2) /= 0 ) &
            &error stop 'center_crop_real3d requires an even nonnegative size difference'
        offset = (padded_box-native_box)/2
        allocate(native(native_box,native_box,native_box), &
            &source=padded(offset+1:offset+native_box, offset+1:offset+native_box, &
            &offset+1:offset+native_box))
    end function center_crop_real3d

    !> Extract one full redundant native Fourier disk from an FFT image.
    function extract_native_fourier_plane( img2d, lims2, sqlp ) result( plane )
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
    end function extract_native_fourier_plane

    !> Gather one KB-weighted value from a packed half-complex Fourier array.
    !! wrap_lower preserves the physical lower bound of the periodic wrap table.
    pure function gather_packed_window( cmat, wrap_lower, wrap, i0, w ) result( value )
        complex, intent(in) :: cmat(:,:,:)
        integer, intent(in) :: wrap_lower
        integer, intent(in) :: wrap(wrap_lower:)
        integer, intent(in) :: i0(3)
        real,    intent(in) :: w(:,:,:)
        complex :: value
        integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
        ny    = size(cmat,2)
        nz    = size(cmat,3)
        value = cmplx(0.,0.)
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
                        value = value+w(di,dj,dk)*cmat(ph,pk,pm)
                    else
                        ph = -hh+1
                        pk = -kk+1; if( -kk < 0 ) pk = pk+ny
                        pm = -mm+1; if( -mm < 0 ) pm = pm+nz
                        value = value+w(di,dj,dk)*conjg(cmat(ph,pk,pm))
                    endif
                end do
            end do
        end do
    end function gather_packed_window

    !> Gather one KB-weighted value and three fixed-cell derivatives in one
    !! packed/Friedel traversal.
    pure subroutine gather_packed_window_grad( cmat, wrap_lower, wrap, i0, w, dw, value, dvalue_dloc )
        complex, intent(in)  :: cmat(:,:,:)
        integer, intent(in)  :: wrap_lower
        integer, intent(in)  :: wrap(wrap_lower:)
        integer, intent(in)  :: i0(3)
        real,    intent(in)  :: w(:,:,:), dw(:,:,:,:)
        complex, intent(out) :: value, dvalue_dloc(3)
        complex :: fcomp
        integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz
        ny          = size(cmat,2)
        nz          = size(cmat,3)
        value       = cmplx(0.,0.)
        dvalue_dloc = cmplx(0.,0.)
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
                    dvalue_dloc = dvalue_dloc+dw(di,dj,dk,:)*fcomp
                end do
            end do
        end do
    end subroutine gather_packed_window_grad

    !> Gather the standard three-tap KB value and fixed-cell derivatives without
    !! materialising three-dimensional weight arrays.
    pure subroutine gather_packed_kb3_window_grad( cmat, wrap_lower, wrap, kbwin, loc, value, dvalue_dloc, inside )
        complex,            intent(in)  :: cmat(:,:,:)
        integer,            intent(in)  :: wrap_lower
        integer,            intent(in)  :: wrap(wrap_lower:)
        type(kbinterpol),   intent(in)  :: kbwin
        real,               intent(in)  :: loc(3)
        complex,            intent(out) :: value, dvalue_dloc(3)
        logical,            intent(out) :: inside
        complex :: fcomp
        real    :: base(3), draw(3), raw(3)
        real    :: dx(3), dy(3), dz(3), wx(3), wy(3), wz(3)
        real    :: sx, sy, sz, dsx, dsy, dsz, wyz, dyz, wydz
        real    :: weight, dweight(3)
        integer :: di, dj, dk, hh, kk, mm, ph, pk, pm, ny, nz, i0(3)

        ! The nearest grid point owns this fixed three-tap cell; derivatives do
        ! not include the discontinuous change to a neighbouring cell.
        i0     = nint(loc) - 1
        inside = .not.(any(i0 < lbound(wrap,1)) .or. any(i0 + 2 > ubound(wrap,1)))
        value       = cmplx(0.,0.)
        dvalue_dloc = cmplx(0.,0.)
        if( .not. inside ) return

        base = real(i0) - loc
        do di = 1, 3
            call kbwin%apod_fast_value_deriv(base + real(di-1), raw, draw)
            wx(di) = raw(1); wy(di) = raw(2); wz(di) = raw(3)
            dx(di) = -draw(1); dy(di) = -draw(2); dz(di) = -draw(3)
        enddo
        ! Normalize each axis separately. The quotient-rule derivatives keep
        ! every 1-D weight sum at one and every derivative sum at zero.
        sx  = sum(wx); sy  = sum(wy); sz  = sum(wz)
        dsx = sum(dx); dsy = sum(dy); dsz = sum(dz)
        wx  = wx / sx; wy = wy / sy; wz = wz / sz
        dx  = (dx - wx * dsx) / sx
        dy  = (dy - wy * dsy) / sy
        dz  = (dz - wz * dsz) / sz

        ! Match the generic oracle's packed/Friedel traversal and accumulation
        ! order while forming separable weights only when they are consumed.
        ny = size(cmat,2)
        nz = size(cmat,3)
        do dk = 1, 3
            mm = wrap(i0(3)+dk-1)
            do dj = 1, 3
                kk   = wrap(i0(2)+dj-1)
                wyz  = wy(dj) * wz(dk)
                dyz  = dy(dj) * wz(dk)
                wydz = wy(dj) * dz(dk)
                do di = 1, 3
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
                    weight     = wx(di) * wyz
                    dweight(1) = dx(di) * wyz
                    dweight(2) = wx(di) * dyz
                    dweight(3) = wx(di) * wydz
                    value       = value + weight * fcomp
                    dvalue_dloc = dvalue_dloc + dweight * fcomp
                enddo
            enddo
        enddo
    end subroutine gather_packed_kb3_window_grad

end module simple_cartesian_fourier
