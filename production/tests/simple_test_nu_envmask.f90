program simple_test_nu_envmask
! Synthetic check of NU-evidence-driven envelope masking.
!
! Ground truth is a centered sphere of band-limited signal shared by both half
! maps, embedded in a larger spherical support that is otherwise independent
! noise. Cross-half consistency therefore exists inside the sphere and nowhere
! else, which is exactly the contrast the evidence margin is supposed to detect.
use simple_core_module_api
use simple_image,     only: image
use simple_image_msk, only: image_msk
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, cleanup_nu_filter, &
    &nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, calc_nu_evidence_margin, &
    &print_nu_envmask_stats, nu_evidence_state, nu_evidence_summary, build_nu_evidence_state, &
    &unpack_nu_evidence_state, get_nu_evidence_summary, nu_evidence_state_is_valid, &
    &print_nu_evidence_summary, assert_nu_evidence_replay_ready, set_nu_solvent_envelope, &
    &NU_EVIDENCE_NBANDS
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"
integer, parameter :: BOX          = 48
real,    parameter :: SMPD         = 2.0
real,    parameter :: MOL_RAD_PX   = 15.0   ! ground-truth molecule radius
real,    parameter :: SUPP_RAD_PX  = 22.0   ! NU support, deliberately generous
real,    parameter :: HARD_RAD_PX  = 18.0   ! emulates a conservative PCG solve support
real,    parameter :: NOISE_SDEV   = 1.0
real,    parameter :: LP_SMOOTH    = 4.0    ! envelope scale; kept small so the
                                            ! assertions test the statistic, not the blur
type(image)             :: even, odd, vol_supp
type(image_msk)         :: envelope
type(nu_envmask_params) :: envp
type(nu_envmask_stats)  :: envstats, envstats_rel
type(nu_evidence_state) :: evstate
type(nu_evidence_summary) :: evsummary
real,    allocatable    :: margin(:), ev_cutoff(:), ev_uncertainty(:), ev_band_support(:,:)
integer, allocatable    :: ev_label(:)
logical, allocatable    :: l_supp(:,:,:), l_env(:,:,:), l_env_rel(:,:,:), l_true(:,:,:), l_hard(:,:,:)
real    :: rmat_even(BOX,BOX,BOX), rmat_odd(BOX,BOX,BOX), sig
real    :: cen, dsq, mean_in, mean_out, recall, fpr, recall_rel, fpr_rel, mskval
real    :: support_in, support_out, null_in, null_out
real    :: unconstrained_null_fraction
integer :: i, j, k, imask, n_true, n_sol, n_hit, n_false, n_ccs, n_ccs_kept, n_clamped
integer :: n_in, n_out
integer(kind=8) :: s1, s2

! ---- synthetic half maps -------------------------------------------------
cen = real(BOX)/2. + 1.
s1  = 1234567_8
s2  = 7654321_8
allocate(l_true(BOX,BOX,BOX), source=.false.)
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            dsq = (real(i)-cen)**2 + (real(j)-cen)**2 + (real(k)-cen)**2
            ! Band-limited common signal, period of roughly five pixels, so an
            ! intermediate bank member wins inside the molecule rather than the
            ! coarsest or the finest one.
            sig = 0.
            if( dsq <= MOL_RAD_PX**2 )then
                l_true(i,j,k) = .true.
                sig = 3.0 * sin(1.26*real(i) + 0.63*real(j)) + &
                     &2.0 * cos(0.94*real(j) - 1.10*real(k)) + &
                     &2.0 * sin(1.05*real(k) + 0.42*real(i))
            endif
            rmat_even(i,j,k) = sig + NOISE_SDEV * lcg_uniform(s1)
            rmat_odd (i,j,k) = sig + NOISE_SDEV * lcg_uniform(s2)
        enddo
    enddo
enddo
call even%new([BOX,BOX,BOX], SMPD)
call odd %new([BOX,BOX,BOX], SMPD)
call even%set_rmat(rmat_even, .false.)
call odd %set_rmat(rmat_odd,  .false.)
call vol_supp%disc([BOX,BOX,BOX], SMPD, SUPP_RAD_PX, l_supp)
if( count(l_supp) < 1 ) THROW_HARD('empty synthetic support mask')
if( .not.all(l_supp .or. .not.l_true) ) THROW_HARD('ground-truth molecule is not inside the support')

! ---- NU setup and evidence ----------------------------------------------
call setup_nu_dmats(even, odd, 2. * SUPP_RAD_PX * SMPD, [real ::], evidence_source='base_unfil')
call optimize_nu_cutoff_finds()
call build_nu_evidence_state(even, odd, evstate)
if( .not.nu_evidence_state_is_valid(evstate) ) THROW_HARD('compact NU evidence state is invalid')
! the PCG replay's readiness contract must hold on this solvent-dominated
! fixture: an inadequate null population here is the zero-null calibration
! failure and must hard-error long before any reconstruction harness runs
call assert_nu_evidence_replay_ready(evstate)
call get_nu_evidence_summary(evstate, evsummary)
unconstrained_null_fraction = evsummary%null_fraction
call unpack_nu_evidence_state(evstate, ev_label, ev_cutoff, ev_uncertainty, ev_band_support)
call print_nu_evidence_summary(evstate)
! The abinitio3D PCG route pins nu_refine=no. Its static evidence state is a
! compatibility boundary: eight signal candidates plus the explicit null,
! four fixed bands, and no adaptive candidate geometry.
if( evsummary%n_bands /= NU_EVIDENCE_NBANDS ) THROW_HARD('static NU evidence band count changed')
if( index(evsummary%provenance, 'adaptive_geometry=') > 0 ) &
    &THROW_HARD('static NU evidence unexpectedly enabled adaptive candidate geometry')
if( evsummary%n_support /= count(l_supp) ) THROW_HARD('compact evidence support size mismatch')
if( evsummary%n_candidates /= 9 ) THROW_HARD('compact evidence must contain null plus the eight-label static bank')
if( abs(evsummary%mskdiam - 2. * SUPP_RAD_PX * SMPD) > TINY ) &
    &THROW_HARD('compact evidence did not retain spherical-support geometry')
if( trim(evsummary%source) /= 'base_unfil' ) THROW_HARD('compact evidence source is not base_unfil')
if( len_trim(evsummary%identity) /= 16 ) THROW_HARD('compact evidence identity is not an FNV-1a hash')
if( any(shape(ev_band_support) /= [count(l_supp),NU_EVIDENCE_NBANDS]) ) &
    &THROW_HARD('compact evidence band-support shape mismatch')
if( any(ev_label < 0) .or. any(ev_label >= evsummary%n_candidates) ) &
    &THROW_HARD('compact evidence selected label is out of range')
if( any(ev_cutoff < 0.) ) THROW_HARD('compact evidence cutoff must be nonnegative')
if( any(ev_uncertainty < 0.) .or. any(ev_uncertainty > 1.) ) &
    &THROW_HARD('compact evidence uncertainty must be in [0,1]')
do i = 2, NU_EVIDENCE_NBANDS
    if( any(ev_band_support(:,i) > ev_band_support(:,i-1) + 1.e-6) ) &
        &THROW_HARD('compact band-support confidence is not monotone coarse-to-fine')
enddo

! The explicit null must win preferentially in independent-noise solvent, while
! coarse-band support must be higher in the common-signal molecule.  These are
! direct gates on the new competitor rather than on the retired binary mask.
support_in  = 0.
support_out = 0.
null_in     = 0.
null_out    = 0.
imask = 0
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            if( .not.l_supp(i,j,k) ) cycle
            imask = imask + 1
            if( l_true(i,j,k) )then
                support_in = support_in + ev_band_support(imask,1)
                if( ev_label(imask) == 0 ) null_in = null_in + 1.
            else
                support_out = support_out + ev_band_support(imask,1)
                if( ev_label(imask) == 0 ) null_out = null_out + 1.
            endif
        enddo
    enddo
enddo
support_in  = support_in / real(count(l_supp .and. l_true))
support_out = support_out / real(count(l_supp .and. .not.l_true))
null_in     = null_in / real(count(l_supp .and. l_true))
null_out    = null_out / real(count(l_supp .and. .not.l_true))
write(logfhandle,'(A,F7.3,A,F7.3)') 'coarse-band support, molecule: ', support_in, '  solvent: ', support_out
write(logfhandle,'(A,F7.3,A,F7.3)') 'null selections, molecule: ', null_in, '  solvent: ', null_out
if( support_in <= support_out ) THROW_HARD('compact NU evidence does not favor signal support inside the molecule')
if( null_out <= null_in ) THROW_HARD('explicit NU null does not preferentially select solvent')
call calc_nu_evidence_margin(margin, LP_SMOOTH, .false.)
if( size(margin) /= count(l_supp) ) THROW_HARD('evidence margin is not packed over the support')
if( any(.not.ieee_is_finite(margin)) ) THROW_HARD('evidence margin produced a non-finite value')
if( any(margin < 0.) ) THROW_HARD('evidence margin must be non-negative by construction')

! The statistic itself must separate before any segmentation is attempted.
mean_in  = 0.
mean_out = 0.
n_in     = 0
n_out    = 0
imask    = 0
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            if( .not.l_supp(i,j,k) ) cycle
            imask = imask + 1
            if( l_true(i,j,k) )then
                mean_in = mean_in + margin(imask)
                n_in    = n_in + 1
            else
                mean_out = mean_out + margin(imask)
                n_out    = n_out + 1
            endif
        enddo
    enddo
enddo
if( n_in < 1 .or. n_out < 1 ) THROW_HARD('degenerate synthetic geometry')
mean_in  = mean_in  / real(n_in)
mean_out = mean_out / real(n_out)
write(logfhandle,'(A,ES12.4,A,ES12.4)') 'mean evidence margin, molecule: ', mean_in, '  solvent: ', mean_out
if( mean_in <= 2. * mean_out ) THROW_HARD('evidence margin does not separate molecule from solvent')

! ---- segmentation --------------------------------------------------------
envp%nsigma      = 3.0
envp%beta        = 1.0
envp%dens_weight = 0.0
envp%lp_smooth   = LP_SMOOTH
envp%l_relative  = .false.
envp%maxits      = 6
call nu_evidence_envelope(envp, l_env, envstats)
call print_nu_envmask_stats(envstats)
if( envstats%n_support /= count(l_supp) ) THROW_HARD('reported support size mismatch')
if( envstats%n_signal < 1 ) THROW_HARD('segmentation produced an empty envelope')
if( any(l_env .and. .not.l_supp) ) THROW_HARD('envelope leaked outside the NU support')

n_true  = count(l_true)
n_sol   = count(l_supp .and. .not.l_true)
n_hit   = count(l_env  .and. l_true)
n_false = count(l_env  .and. l_supp .and. .not.l_true)
recall  = real(n_hit)   / real(n_true)
fpr     = real(n_false) / real(n_sol)
write(logfhandle,'(A,F7.3,A,F7.3)') 'envelope recall: ', recall, '  solvent false-positive rate: ', fpr
if( recall < 0.70 ) THROW_HARD('envelope recovered too little of the ground-truth molecule')
if( fpr    > 0.20 ) THROW_HARD('envelope admitted too much solvent')

! ---- the scale-free margin must recover the same molecule ----------------
call calc_nu_evidence_margin(margin, LP_SMOOTH, .true.)
if( any(.not.ieee_is_finite(margin)) ) THROW_HARD('scale-free evidence produced a non-finite value')
if( any(margin < 0.) ) THROW_HARD('scale-free evidence must be non-negative by construction')
mean_in  = 0.
mean_out = 0.
imask    = 0
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            if( .not.l_supp(i,j,k) ) cycle
            imask = imask + 1
            if( l_true(i,j,k) )then
                mean_in = mean_in + margin(imask)
            else
                mean_out = mean_out + margin(imask)
            endif
        enddo
    enddo
enddo
mean_in  = mean_in  / real(n_in)
mean_out = mean_out / real(n_out)
write(logfhandle,'(A,ES12.4,A,ES12.4,A,ES12.4,A,ES12.4)') &
    &'mean relative margin, molecule: ', mean_in, '  solvent: ', mean_out, &
    &'  range: ', minval(margin), ' to ', maxval(margin)
if( mean_in <= 2. * mean_out ) THROW_HARD('scale-free evidence does not separate molecule from solvent')
envp%l_relative = .true.
call nu_evidence_envelope(envp, l_env_rel, envstats_rel)
call print_nu_envmask_stats(envstats_rel)
if( envstats_rel%n_signal < 1 ) THROW_HARD('scale-free segmentation produced an empty envelope')
if( any(l_env_rel .and. .not.l_supp) ) THROW_HARD('scale-free envelope leaked outside the NU support')
recall_rel = real(count(l_env_rel .and. l_true))                   / real(n_true)
fpr_rel    = real(count(l_env_rel .and. l_supp .and. .not.l_true)) / real(n_sol)
write(logfhandle,'(A,F7.3,A,F7.3)') 'scale-free recall: ', recall_rel, '  solvent false-positive rate: ', fpr_rel
if( recall_rel < 0.70 ) THROW_HARD('scale-free envelope recovered too little of the ground-truth molecule')
if( fpr_rel    > 0.20 ) THROW_HARD('scale-free envelope admitted too much solvent')
envp%l_relative = .false.

! ---- topology and morphology tail ---------------------------------------
call envelope%envmask3D_from_lmask(l_env, SMPD, 1, 3, 0.1, .true., n_ccs, n_ccs_kept)
if( n_ccs      < 1 ) THROW_HARD('no connected components in the segmented envelope')
if( n_ccs_kept < 1 ) THROW_HARD('component size filtering discarded every component')
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            mskval = envelope%get([i,j,k])
            if( .not.ieee_is_finite(mskval) ) THROW_HARD('soft envelope mask has a non-finite value')
            if( mskval < -1.e-5 .or. mskval > 1.+1.e-5 ) THROW_HARD('soft envelope mask is not in [0,1]')
        enddo
    enddo
enddo
if( envelope%get([nint(cen),nint(cen),nint(cen)]) < 0.5 ) THROW_HARD('molecule centre is not inside the soft envelope')
if( envelope%get([1,1,1]) > 1.e-5 ) THROW_HARD('box corner is inside the soft envelope')

! ---- PCG evidence background constraint ---------------------------------
! The readiness null fraction remains that of the unconstrained spherical
! Potts field, while the replay state fixes envelope-background voxels to the
! coarsest bank member: coarse support 1, every finer support 0.
call set_nu_solvent_envelope(envelope, source='test_density_envelope')
call build_nu_evidence_state(even, odd, evstate)
call assert_nu_evidence_replay_ready(evstate)
call get_nu_evidence_summary(evstate, evsummary)
if( abs(evsummary%null_fraction - unconstrained_null_fraction) > 1.e-6 ) &
    &THROW_HARD('envelope constraint changed the broad-support NU null readiness statistic')
deallocate(ev_label, ev_cutoff, ev_uncertainty, ev_band_support)
call unpack_nu_evidence_state(evstate, ev_label, ev_cutoff, ev_uncertainty, ev_band_support)
n_clamped = 0
imask = 0
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            if( .not.l_supp(i,j,k) ) cycle
            imask = imask + 1
            if( envelope%get([i,j,k]) >= 0.5 ) cycle
            n_clamped = n_clamped + 1
            if( ev_label(imask) /= 1 ) THROW_HARD('NU background was not assigned the coarsest signal label')
            if( ev_cutoff(imask) <= TINY ) THROW_HARD('NU background coarsest assignment has no low-pass cutoff')
            if( abs(ev_band_support(imask,1) - 1.) > 1.e-6 ) &
                &THROW_HARD('NU background does not preserve the coarsest evidence band')
            if( maxval(abs(ev_band_support(imask,2:NU_EVIDENCE_NBANDS))) > 1.e-6 ) &
                &THROW_HARD('NU background incorrectly supports detail finer than the coarsest band')
        enddo
    enddo
enddo
if( n_clamped < 1 ) THROW_HARD('synthetic NU envelope produced no constrained background voxels')

call cleanup_nu_filter()
if( .not.nu_evidence_state_is_valid(evstate) ) &
    &THROW_HARD('compact NU evidence changed when mutable unary state was released')

! A density-supported PCG pair is identically zero outside the conservative
! solve envelope while the NU evidence domain remains the broader sphere. The
! exact-zero boundary must not collapse the radial noise estimate or overflow
! the explicit zero-predictor unary during evidence compaction.
allocate(l_hard(BOX,BOX,BOX), source=.false.)
do k = 1,BOX
    do j = 1,BOX
        do i = 1,BOX
            dsq = (real(i)-cen)**2 + (real(j)-cen)**2 + (real(k)-cen)**2
            l_hard(i,j,k) = dsq <= HARD_RAD_PX**2
            if( .not.l_hard(i,j,k) )then
                rmat_even(i,j,k) = 0.
                rmat_odd (i,j,k) = 0.
            endif
        enddo
    enddo
enddo
if( count(l_supp .and. .not.l_hard) < 1 ) THROW_HARD('hard-support regression has no exact-zero boundary')
call even%set_rmat(rmat_even, .false.)
call odd %set_rmat(rmat_odd,  .false.)
call setup_nu_dmats(even, odd, 2. * SUPP_RAD_PX * SMPD, [real ::], evidence_source='base_unfil')
call optimize_nu_cutoff_finds()
call build_nu_evidence_state(even, odd, evstate)
if( .not.nu_evidence_state_is_valid(evstate) ) &
    &THROW_HARD('hard-supported pair produced invalid NU evidence')
call assert_nu_evidence_replay_ready(evstate)
call cleanup_nu_filter()

call even%kill
call odd%kill
call vol_supp%kill
call envelope%kill_bimg
deallocate(margin, ev_label, ev_cutoff, ev_uncertainty, ev_band_support, l_supp, l_env, l_true, l_hard)
if( allocated(l_env_rel) ) deallocate(l_env_rel)
write(logfhandle,'(A)') 'NU evidence envelope masking test passed'

contains

    !>  Deterministic uniform deviate on [-0.5,0.5]. Two separately seeded streams
    !!  give the half maps independent noise without depending on the platform RNG.
    real function lcg_uniform( seed )
        integer(kind=8), intent(inout) :: seed
        seed = mod(1103515245_8 * seed + 12345_8, 2147483648_8)
        lcg_uniform = real(seed) / 2147483648.0 - 0.5
    end function lcg_uniform

end program simple_test_nu_envmask
