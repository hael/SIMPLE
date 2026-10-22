!@descr: neutral Cartesian Fourier lattice embedding, extraction, and packed gathers
module simple_cartesian_fourier
use simple_defs,  only: PI
use simple_image, only: image
implicit none
private

public :: center_embed_real3d, center_crop_real3d
public :: extract_native_fourier_plane
public :: gather_packed_window, gather_packed_window_grad, gather_packed_bspline3_window_grad
public :: bspline3_stencil_envelope_1d

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

    !> Gather one value and its three derivatives with the separable quadratic B-spline on the
    !! three taps around the nearest grid point. Its weights and slopes vanish at the taps that
    !! enter or leave as the nearest point changes, so value and derivatives are continuous.
    pure subroutine gather_packed_bspline3_window_grad( cmat, wrap_lower, wrap, loc, value, dvalue_dloc, inside )
        complex, intent(in)  :: cmat(:,:,:)
        integer, intent(in)  :: wrap_lower
        integer, intent(in)  :: wrap(wrap_lower:)
        real,    intent(in)  :: loc(3)
        complex, intent(out) :: value, dvalue_dloc(3)
        logical, intent(out) :: inside
        complex :: f(3), fx, fdx, fy(3), fdy(3), fxdy(3)
        real    :: t(3), wx(3), wy(3), wz(3), dx(3), dy(3), dz(3)
        integer :: di, dj, dk, hh(3), kk(3), mm(3), ph(3), pk(3), pm(3), pkc(3), pmc(3), ny, nz, i0(3)
        logical :: l_direct(3), l_run
        i0     = nint(loc) - 1
        inside = .not.(any(i0 < lbound(wrap,1)) .or. any(i0 + 2 > ubound(wrap,1)))
        value       = cmplx(0.,0.)
        dvalue_dloc = cmplx(0.,0.)
        if( .not. inside ) return
        ! offset t in [-1/2,1/2] from the nearest point; taps at distances -1-t, -t and 1-t, whose
        ! weights sum to one and whose derivatives in loc sum to zero for every t
        t  = loc - real(i0 + 1)
        wx = [0.5*(t(1)-0.5)**2, 0.75-t(1)**2, 0.5*(t(1)+0.5)**2]
        wy = [0.5*(t(2)-0.5)**2, 0.75-t(2)**2, 0.5*(t(2)+0.5)**2]
        wz = [0.5*(t(3)-0.5)**2, 0.75-t(3)**2, 0.5*(t(3)+0.5)**2]
        dx = [t(1)-0.5, -2.*t(1), t(1)+0.5]
        dy = [t(2)-0.5, -2.*t(2), t(2)+0.5]
        dz = [t(3)-0.5, -2.*t(3), t(3)+0.5]
        ! packed/Friedel indices per axis: a tap with h < 0 reads the conjugate at (-h,-k,-m)
        ny = size(cmat,2)
        nz = size(cmat,3)
        do di = 1, 3
            hh(di) = wrap(i0(1)+di-1)
            kk(di) = wrap(i0(2)+di-1)
            mm(di) = wrap(i0(3)+di-1)
        end do
        l_direct = hh >= 0
        l_run    = all(l_direct) .and. hh(3) == hh(1) + 2
        ph  = abs(hh) + 1
        pk  =  kk + 1; where(  kk < 0 ) pk  = pk  + ny
        pm  =  mm + 1; where(  mm < 0 ) pm  = pm  + nz
        pkc = -kk + 1; where( -kk < 0 ) pkc = pkc + ny
        pmc = -mm + 1; where( -mm < 0 ) pmc = pmc + nz
        ! separable contraction: x, then y, then z
        do dk = 1, 3
            fy(dk)   = cmplx(0.,0.)
            fdy(dk)  = cmplx(0.,0.)
            fxdy(dk) = cmplx(0.,0.)
            do dj = 1, 3
                if( l_run )then
                    f = cmat(ph(1):ph(3),pk(dj),pm(dk))
                else
                    do di = 1, 3
                        if( l_direct(di) )then
                            f(di) = cmat(ph(di),pk(dj),pm(dk))
                        else
                            f(di) = conjg(cmat(ph(di),pkc(dj),pmc(dk)))
                        endif
                    end do
                endif
                fx  = wx(1)*f(1) + wx(2)*f(2) + wx(3)*f(3)
                fdx = dx(1)*f(1) + dx(2)*f(2) + dx(3)*f(3)
                fy(dk)   = fy(dk)   + wy(dj)*fx
                fdy(dk)  = fdy(dk)  + dy(dj)*fx
                fxdy(dk) = fxdy(dk) + wy(dj)*fdx
            end do
        end do
        value          = wz(1)*fy(1)   + wz(2)*fy(2)   + wz(3)*fy(3)
        dvalue_dloc(1) = wz(1)*fxdy(1) + wz(2)*fxdy(2) + wz(3)*fxdy(3)
        dvalue_dloc(2) = wz(1)*fdy(1)  + wz(2)*fdy(2)  + wz(3)*fdy(3)
        dvalue_dloc(3) = dz(1)*fy(1)   + dz(2)*fy(2)   + dz(3)*fy(3)
    end subroutine gather_packed_bspline3_window_grad

    !> Real-space response along one axis of a box of n of the quadratic B-spline stencil at
    !! on-grid samples, weights 1/8, 3/4, 1/8: 3/4 + cos(2 pi (x-c)/n)/4, centre c = n/2+1.
    subroutine bspline3_stencil_envelope_1d( n, env1d )
        integer,           intent(in)  :: n
        real, allocatable, intent(out) :: env1d(:)
        integer :: x
        allocate(env1d(n))
        do x = 1, n
            env1d(x) = 0.75 + 0.25*cos(2.*PI*real(x - (n/2 + 1))/real(n))
        end do
    end subroutine bspline3_stencil_envelope_1d

end module simple_cartesian_fourier
