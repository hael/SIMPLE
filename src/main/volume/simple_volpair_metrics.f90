!@descr: two maps of one frame and sampling compared inside a soft mask: their FSC and their correlation up to a resolution
! compare_volpair reads both maps, soft-masks them (mask diameter in A, zero
! background) and returns the FSC between them with its resolution axis, and the
! correlation of their Fourier coefficients over a band (hp..lp A), which is the
! real-space correlation of the band-passed masked maps. The maps must agree in
! box and sampling; ok is false when they do not, when a file is missing, or
! when a value is not finite.
module simple_volpair_metrics
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_image,   only: image
use simple_imghead, only: find_ldim_nptcls, find_img_smpd
implicit none

public :: compare_volpair
private
#include "simple_local_flags.inc"

real, parameter :: SMPD_TOL = 1.e-3 !< sampling agreement required of the two maps (A)

contains

    !> FSC and band correlation of the maps in fname1 and fname2, both
    !! soft-masked at mskdiam (A); the band runs from hp (A, default the second
    !! Fourier shell) to lp (A; lp <= 0 means Nyquist)
    subroutine compare_volpair( fname1, fname2, mskdiam, lp, corr, fsc, res, ok, hp )
        class(string),     intent(in)  :: fname1, fname2
        real,              intent(in)  :: mskdiam, lp
        real,              intent(out) :: corr
        real, allocatable, intent(out) :: fsc(:), res(:)
        logical,           intent(out) :: ok
        real, optional,    intent(in)  :: hp
        type(image) :: vol1, vol2
        real    :: smpd1, smpd2
        integer :: ldim1(3), ldim2(3), nsections
        corr = 0.
        ok   = .false.
        if( .not. file_exists(fname1) .or. .not. file_exists(fname2) ) return
        call find_ldim_nptcls(fname1, ldim1, nsections)
        call find_ldim_nptcls(fname2, ldim2, nsections)
        smpd1 = find_img_smpd(fname1)
        smpd2 = find_img_smpd(fname2)
        if( any(ldim1 /= ldim2) .or. abs(smpd1 - smpd2) > SMPD_TOL ) return
        call vol1%new(ldim1, smpd1)
        call vol2%new(ldim1, smpd1)
        call vol1%read(fname1)
        call vol2%read(fname2)
        call vol1%mask3D_soft(0.5 * mskdiam / smpd1, backgr=0.)
        call vol2%mask3D_soft(0.5 * mskdiam / smpd1, backgr=0.)
        call vol1%fft()
        call vol2%fft()
        allocate(fsc(vol1%get_filtsz()), source=0.)
        call vol1%fsc(vol2, fsc)
        res = vol1%get_res()
        if( lp > 0. )then
            if( present(hp) )then
                corr = vol1%corr(vol2, lp_dyn=lp, hp_dyn=hp)
            else
                corr = vol1%corr(vol2, lp_dyn=lp)
            endif
        else
            if( present(hp) )then
                corr = vol1%corr(vol2, hp_dyn=hp)
            else
                corr = vol1%corr(vol2)
            endif
        endif
        ok = all(ieee_is_finite(fsc)) .and. ieee_is_finite(corr)
        call vol1%kill
        call vol2%kill
    end subroutine compare_volpair

end module simple_volpair_metrics
