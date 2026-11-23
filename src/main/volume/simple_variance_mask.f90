!@descr: soft mask from the population-weighted cross-half variance of state half maps
module simple_variance_mask
use simple_core_module_api, only: dp, find_img_smpd, find_ldim_nptcls, int2str, logfhandle, PI, string
use simple_image,           only: image
implicit none

public :: write_variance_mask
private

real, parameter :: VARMASK_FRAC = 0.2, VARMASK_SUPPORT = 0.1, VARMASK_DIL = 3., VARMASK_FALL = 6., VARMASK_LP = 15., VARMASK_SIG = 1.5

contains

    subroutine write_variance_mask( evens, odds, w, mskdiam, fname )
        type(string),     intent(in) :: evens(:), odds(:)
        real(dp),         intent(in) :: w(:)
        real,             intent(in) :: mskdiam
        character(len=*), intent(in) :: fname
        type(image)          :: img
        real,    allocatable :: vs(:,:,:), m(:,:,:)
        logical, allocatable :: sph(:,:,:), sup(:,:,:), core(:,:,:)
        integer :: n
        real    :: smpd
        call half_variance(evens, odds, w, mskdiam, n, smpd, sph, sup, vs)
        core = vs > VARMASK_FRAC * maxval(vs)
        if( .not. any(core) )then
            m = merge(1., 0., sph)
        else
            m    = dist_vox(core, n, smpd)
            core = m <= VARMASK_DIL
            m    = dist_vox(core, n, smpd)
            m    = exp(-m**2 / (2. * VARMASK_FALL**2))
            where( .not. sph ) m = 0.
        endif
        call img%new([n,n,n], smpd)
        call img%set_rmat(m, .false.)
        call img%write(string(fname))
        call img%set_rmat(vs, .false.)
        call img%write(string(fname(1:len(fname)-4)//'_variance.mrc'))
        call img%kill
        write(logfhandle,'(A)') '>>> VARIANCE MASK '//fname//' from '//int2str(size(evens))//' state maps: '//&
            &int2str(count(m > 0.5))//' voxels above 0.5 ('//int2str(nint(100. * real(count(m > 0.5 .and. sup)) / &
            &real(max(1,count(sup)))))//' % of the molecule support)'
    end subroutine write_variance_mask

    subroutine half_variance( evens, odds, w, mskdiam, n, smpd, sph, sup, vs )
        type(string),         intent(in)  :: evens(:), odds(:)
        real(dp),             intent(in)  :: w(:)
        real,                 intent(in)  :: mskdiam
        integer,              intent(out) :: n
        real,                 intent(out) :: smpd
        logical, allocatable, intent(out) :: sph(:,:,:), sup(:,:,:)
        real,    allocatable, intent(out) :: vs(:,:,:)
        real,     allocatable :: e(:,:,:,:), o(:,:,:,:), em(:,:,:), om(:,:,:), v(:,:,:)
        real(dp), allocatable :: wn(:)
        integer :: ldim(3), nm, j, ifoo
        real    :: sc, fwhm
        call find_ldim_nptcls(evens(1), ldim, ifoo)
        n    = ldim(1)
        nm   = size(evens)
        smpd = find_img_smpd(evens(1))
        wn   = w / sum(w)
        fwhm = VARMASK_LP * 2. * log(2.) / PI
        allocate(e(n,n,n,nm), o(n,n,n,nm))
        do j = 1, nm
            e(:,:,:,j) = read_smoothed(evens(j), n, smpd, fwhm)
            o(:,:,:,j) = read_smoothed(odds(j),  n, smpd, fwhm)
        end do
        sph = sphere(n, mskdiam / 2. / smpd)
        em  = wsum(e)
        om  = wsum(o)
        sup = (0.5 * (em + om) > VARMASK_SUPPORT * maxval(0.5 * (em + om), mask=sph)) .and. sph
        do j = 1, nm
            sc = sum(e(:,:,:,j) * em, mask=sup) / max(sum(e(:,:,:,j)**2, mask=sup), 1.e-30)
            e(:,:,:,j) = e(:,:,:,j) * sc
            sc = sum(o(:,:,:,j) * om, mask=sup) / max(sum(o(:,:,:,j)**2, mask=sup), 1.e-30)
            o(:,:,:,j) = o(:,:,:,j) * sc
        end do
        em = wsum(e)
        om = wsum(o)
        allocate(v(n,n,n), source=0.)
        do j = 1, nm
            v = v + real(wn(j)) * (e(:,:,:,j) - em) * (o(:,:,:,j) - om)
        end do
        vs = smoothed(v, n, smpd, 2.35482 * VARMASK_SIG * smpd)
        where( .not. sph ) vs = 0.

    contains

        function wsum( x ) result( s )
            real, intent(in)  :: x(:,:,:,:)
            real, allocatable :: s(:,:,:)
            integer :: jj
            allocate(s(n,n,n), source=0.)
            do jj = 1, nm
                s = s + real(wn(jj)) * x(:,:,:,jj)
            end do
        end function wsum

    end subroutine half_variance

    function read_smoothed( fname, n, smpd, fwhm ) result( x )
        type(string), intent(in) :: fname
        integer,      intent(in) :: n
        real,         intent(in) :: smpd, fwhm
        real, allocatable :: x(:,:,:)
        type(image) :: vol
        call vol%new([n,n,n], smpd)
        call vol%read(fname)
        call vol%fft
        call vol%bpgau3D(0., fwhm)
        call vol%ifft
        x = vol%get_rmat()
        call vol%kill
    end function read_smoothed

    function smoothed( x, n, smpd, fwhm ) result( y )
        real,    intent(in) :: x(:,:,:), smpd, fwhm
        integer, intent(in) :: n
        real, allocatable :: y(:,:,:)
        type(image) :: vol
        call vol%new([n,n,n], smpd)
        call vol%set_rmat(x, .false.)
        call vol%fft
        call vol%bpgau3D(0., fwhm)
        call vol%ifft
        y = vol%get_rmat()
        call vol%kill
    end function smoothed

    function dist_vox( mk, n, smpd ) result( d )
        logical, intent(in) :: mk(:,:,:)
        integer, intent(in) :: n
        real,    intent(in) :: smpd
        real, allocatable :: d(:,:,:)
        type(image) :: vol
        call vol%new([n,n,n], smpd)
        call vol%set_rmat(merge(0., 1., mk), .false.)
        call vol%distance_transform
        d = vol%get_rmat() / smpd
        call vol%kill
    end function dist_vox

    function sphere( n, rad ) result( sph )
        integer, intent(in) :: n
        real,    intent(in) :: rad
        logical, allocatable :: sph(:,:,:)
        integer :: i, j, k, c
        c = n / 2 + 1
        allocate(sph(n,n,n))
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    sph(i,j,k) = real((i-c)**2 + (j-c)**2 + (k-c)**2) <= rad**2
                end do
            end do
        end do
    end function sphere

end module simple_variance_mask
