program simple_test_ctf
use simple_core_module_api
use simple_image,           only: image
use simple_ctf,             only: ctf
use simple_memoize_ft_maps, only: memoize_ft_maps
use simple_test_utils,      only: assert_real, assert_true, report_summary, tests_failed
implicit none

integer, parameter :: LDIM(3) = [64,64,1]
integer, parameter :: H = 5, K = 3
real,    parameter :: SMPD = 1.0, DFX = 2.0, DFY = 2.3, ANGAST = 27.
real,    parameter :: KV = 300., CS = 2.0, AC = 0.1, TOL = 2.e-6
real,    parameter :: PHASES(3) = [0., PI/4., PIO2]
type(image)        :: img
type(ctf)          :: tfun
type(ctfparams)    :: ctfparms
type(ctfvars)      :: ctfvals
type(fplane_type)  :: fplane
real               :: spa_freq_sq, ang, expected, actual, conventional, shifted
integer            :: iph, hp, kp

call img%new(LDIM, SMPD)
tfun = ctf(SMPD, KV, CS, AC)
call memoize_ft_maps(LDIM, SMPD)
spa_freq_sq = (real(H)/real(LDIM(1)))**2 + (real(K)/real(LDIM(2)))**2
ang         = atan2(real(K),real(H))

call assert_real(0., canonical_phshift(PI), TOL, &
    &'phase canonicalization maps pi to zero')
call assert_real(3.*PI/4., canonical_phshift(-PI/4.), TOL, &
    &'phase canonicalization maps negative values into [0,pi)')
call tfun%init(DFX, DFY, ANGAST)
ctfvals = tfun%get_ctfvars(PI + PI/4.)
call assert_real(PI/4., ctfvals%phshift, TOL, &
    &'ctfvars carries the canonical additive phase shift')
call assert_real(tfun%eval(spa_freq_sq, ang, PI/4.), &
    &tfun%eval(spa_freq_sq, ang, PI + PI/4.), TOL, &
    &'scalar CTF evaluation consumes the canonical phase')

! The image-layer fast kernel must agree with the scalar CTF for conventional,
! intermediate, and near-quadrature phase shifts.
do iph = 1,size(PHASES)
    call img%ctf2img(tfun, DFX, DFY, ANGAST, PHASES(iph))
    call tfun%init(DFX, DFY, ANGAST)
    expected = tfun%eval(spa_freq_sq, ang, PHASES(iph))
    actual   = real(img%get_fcomp2D(H,K))
    call assert_real(expected, actual, TOL, 'fast CTF kernel agrees with scalar phase-shift CTF')
    if( iph == 1 ) conventional = actual
    if( iph == size(PHASES) ) shifted = actual
enddo
call img%ctf2img(tfun, DFX, DFY, ANGAST, 0.)
call assert_real(conventional, real(img%get_fcomp2D(H,K)), TOL, &
    &'explicit zero phase shift preserves conventional CTF behavior')
call assert_true(abs(shifted-conventional) > 0.1, &
    &'nonzero phase shift materially changes the transfer function')

! Reconstruction uses H*y and H^2 for raw observations, and |H|*y and H^2
! for phase-flipped observations.  A unit Fourier image makes both contracts
! directly observable at a native-grid sample.
call img%set_cmat(cmplx(1.,0.))
ctfparms%smpd        = SMPD
ctfparms%kv          = KV
ctfparms%cs          = CS
ctfparms%fraca       = AC
ctfparms%dfx         = DFX
ctfparms%dfy         = DFY
ctfparms%angast      = ANGAST
ctfparms%phshift     = PI/4.
ctfparms%ctfflag     = CTFFLAG_YES
call img%gen_fplane4rec([0,LDIM(1)/(2*OSMPL_PAD_FAC)], SMPD, ctfparms, &
    &[0.,0.], fplane)
hp = OSMPL_PAD_FAC*H
kp = -OSMPL_PAD_FAC*K
actual = real(fplane%cmplx_plane(hp,kp))
call assert_real(actual*actual, fplane%ctfsq_plane(hp,kp), 5.e-6, &
    &'raw 3D restoration stores H*y over H^2')

ctfparms%ctfflag = CTFFLAG_FLIP
call img%gen_fplane4rec([0,LDIM(1)/(2*OSMPL_PAD_FAC)], SMPD, ctfparms, &
    &[0.,0.], fplane)
actual = real(fplane%cmplx_plane(hp,kp))
call assert_true(actual >= 0., 'phase-flipped 3D restoration uses nonnegative |H|')
call assert_real(actual*actual, fplane%ctfsq_plane(hp,kp), 5.e-6, &
    &'phase-flipped 3D restoration stores |H|*y over H^2')

call report_summary()
call img%kill
if( tests_failed > 0 ) error stop 1
end program simple_test_ctf
