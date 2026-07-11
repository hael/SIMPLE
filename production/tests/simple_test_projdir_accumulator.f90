program simple_test_projdir_accumulator
use simple_core_module_api
use simple_classaverager, only: fourier_2d_accumulator
use simple_ftiter,       only: ftiter
use simple_test_utils,   only: assert_true
implicit none

integer, parameter :: BOX = 32
real,    parameter :: SIGNAL_SCALE = 2.0
type(fourier_2d_accumulator) :: sums
type(fplane_type) :: padded, compact
type(ftiter) :: fit_pd
integer :: lims_pd(3,2), h, k, hp, kp, sh
real    :: rho, max_num_err, max_imag, rho_sum

fit_pd = ftiter([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], 1.0)
lims_pd = fit_pd%loop_lims(3)
allocate(padded%cmplx_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0), source=cmplx(0.,0.))
allocate(padded%ctfsq_plane(lims_pd(1,1):lims_pd(1,2),lims_pd(2,1):0), source=0.0)
padded%frlims = lims_pd
padded%nyq    = fit_pd%get_lfny(1)

! Populate only the native samples carried by gen_fplane4rec.  A real-valued
! numerator proportional to CTF^2 makes the expected result independent of
! rotation while still exercising the complete 2D KB splat and export path.
do k = -BOX/2,0
    kp = OSMPL_PAD_FAC * k
    do h = -BOX/2,BOX/2
        hp = OSMPL_PAD_FAC * h
        sh = nint(hyp(h,k))
        if( sh > BOX/2 ) cycle
        rho = 0.25 + real(mod(abs(3*h+5*k),11)) / 11.0
        padded%ctfsq_plane(hp,kp) = rho
        padded%cmplx_plane(hp,kp) = cmplx(SIGNAL_SCALE*rho,0.0)
    enddo
enddo

call sums%new([BOX,BOX], 1)
call sums%add_fplane(37.25, padded, 1)
call sums%export_fplane(1, compact)

max_num_err = maxval(abs(real(compact%cmplx_plane) - &
    &(OSMPL_PAD_FAC**2)*SIGNAL_SCALE*compact%ctfsq_plane))
max_imag = maxval(abs(aimag(compact%cmplx_plane)))
rho_sum  = sum(compact%ctfsq_plane)
call assert_true(rho_sum > 0.0, 'compact projection-direction CTF^2 sum is populated')
call assert_true(max_num_err < 2.e-4, 'compact numerator and CTF^2 use identical normalized KB weights')
call assert_true(max_imag < 2.e-5, 'real Friedel-symmetric input remains real after compact export')
call assert_true(ubound(compact%cmplx_plane,2) == 0, 'compact export stores only the k<=0 half-plane')
call assert_true(compact%nyq == BOX/2, 'compact export uses the native-grid Nyquist limit')

call sums%kill
if( allocated(padded%cmplx_plane) ) deallocate(padded%cmplx_plane)
if( allocated(padded%ctfsq_plane) ) deallocate(padded%ctfsq_plane)
if( allocated(compact%cmplx_plane) ) deallocate(compact%cmplx_plane)
if( allocated(compact%ctfsq_plane) ) deallocate(compact%ctfsq_plane)
end program simple_test_projdir_accumulator
