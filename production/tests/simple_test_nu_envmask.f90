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
    &print_nu_envmask_stats
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none
#include "simple_local_flags.inc"
integer, parameter :: BOX          = 48
real,    parameter :: SMPD         = 2.0
real,    parameter :: MOL_RAD_PX   = 15.0   ! ground-truth molecule radius
real,    parameter :: SUPP_RAD_PX  = 22.0   ! NU support, deliberately generous
real,    parameter :: NOISE_SDEV   = 1.0
real,    parameter :: LP_SMOOTH    = 4.0    ! envelope scale; kept small so the
                                            ! assertions test the statistic, not the blur
type(image)             :: even, odd, vol_supp
type(image_msk)         :: envelope
type(nu_envmask_params) :: envp
type(nu_envmask_stats)  :: envstats, envstats_rel
real,    allocatable    :: margin(:)
logical, allocatable    :: l_supp(:,:,:), l_env(:,:,:), l_env_rel(:,:,:), l_true(:,:,:)
real    :: rmat_even(BOX,BOX,BOX), rmat_odd(BOX,BOX,BOX), sig
real    :: cen, dsq, mean_in, mean_out, recall, fpr, recall_rel, fpr_rel, mskval
integer :: i, j, k, imask, n_true, n_sol, n_hit, n_false, n_ccs, n_ccs_kept
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
call setup_nu_dmats(even, odd, 2. * SUPP_RAD_PX * SMPD, [real ::])
call optimize_nu_cutoff_finds()
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

call cleanup_nu_filter()
call even%kill
call odd%kill
call vol_supp%kill
call envelope%kill_bimg
deallocate(margin, l_supp, l_env, l_true)
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
