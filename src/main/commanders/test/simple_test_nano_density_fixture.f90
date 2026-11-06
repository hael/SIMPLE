!@descr: the transfer a rendered nanoparticle was made with, as a perfect measurement returns it: a fixture shared by the nano density tester and the SINGLE high-level tests (simple_nano_density)
module simple_test_nano_density_fixture
use simple_defs,         only: dp
use simple_nano_density, only: render_atoms
implicit none
private
public :: rendered_transfer

contains

    ! the transfer a particle was rendered with, as a perfect measurement returns it: the atoms (xyz in voxels, ints,
    ! bfacs) rendered without a transfer in a box of nbox at smpd, then per Fourier shell of width ds (nsh shells from
    ! the origin, the origin left out) the notch 1 - amp exp(-bh s**2/4) averaged over the shell's coefficients with
    ! their power |F|**2 as weights (t_ref); have marks the shells holding at least 1% of the peak power
    subroutine rendered_transfer( xyz, ints, bfacs, nbox, smpd_r, amp, bh, ds, nsh, t_ref, have )
        use simple_image, only: image
        real,                 intent(in)  :: xyz(:,:), ints(:), bfacs(:), smpd_r, amp, bh, ds
        integer,              intent(in)  :: nbox, nsh
        real,    allocatable, intent(out) :: t_ref(:)
        logical, allocatable, intent(out) :: have(:)
        type(image) :: img
        real,     allocatable :: arr(:,:,:)
        complex,  allocatable :: c(:,:,:)
        real(dp), allocatable :: pw(:), pwt(:)
        integer :: cs(3), i, j, l, h, k, m, ish
        real    :: sf, w
        allocate(arr(nbox,nbox,nbox))
        call render_atoms(arr, smpd_r, xyz, ints, bfacs, 0., 0.)
        call img%new([nbox,nbox,nbox], smpd_r, wthreads=.false.)
        call img%set_rmat(arr, .false.)
        call img%fft()
        c  = img%get_cmat()
        cs = shape(c)
        call img%kill
        allocate(pw(nsh), pwt(nsh), source=0._dp)
        do l = 1,cs(3)
            m = l - 1
            if( m > nbox / 2 ) m = m - nbox
            do j = 1,cs(2)
                k = j - 1
                if( k > nbox / 2 ) k = k - nbox
                do i = 1,cs(1)
                    h = i - 1
                    if( h == 0 .and. k == 0 .and. m == 0 ) cycle
                    sf  = sqrt(real(h**2 + k**2 + m**2)) / (real(nbox) * smpd_r)
                    ish = min(int(sf / ds) + 1, nsh)
                    w   = 2.
                    if( h == 0 .or. 2 * h == nbox ) w = 1.
                    pw(ish)  = pw(ish)  + real(w * abs(c(i,j,l))**2, dp)
                    pwt(ish) = pwt(ish) + real(w * abs(c(i,j,l))**2 * (1. - amp * exp(-bh * sf**2 / 4.)), dp)
                enddo
            enddo
        enddo
        have  = pw >= 0.01_dp * maxval(pw) .and. pw > 0._dp
        t_ref = real(pwt / max(pw, tiny(1._dp)))
    end subroutine rendered_transfer

end module simple_test_nano_density_fixture
