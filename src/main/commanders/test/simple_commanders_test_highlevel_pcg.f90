!@descr: PCG operator and solver gate of simple_commanders_test_highlevel
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_pcg
implicit none
#include "simple_local_flags.inc"

contains

!>  \brief  Fail-fast gate of the PCG operator and solver: stages 1-14 plus 3b; stages after a failure are skipped.
!  Stage table: doc/policies/3D/reconstruct3D_pcg_policy.md sec. 9. Operator data are forward_plane projections
!  (inverse crime: algebra only); stage 8 projects phantom/env and is the only envelope check.
module subroutine exec_test_pcg_recon( self, cline )
    use simple_reconstructor_pcg, only: reconstructor_pcg, pcg_solver_outcome, &
        &PCG_OP_MATRIXFREE, PCG_OP_KERNEL, PCG_STOP_INDEFINITE
    use simple_sym,               only: sym
    use simple_image,             only: image
    use simple_matcher_ptcl_io,   only: prep_rec_observation
    class(commander_test_pcg_recon), intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    integer,          parameter :: BOX = 32, CROP_BOX = 16, NPROJS = 40, NBLOBS = 4, NCTF = 5
    integer,          parameter :: BATCHSZ = 7          ! deliberately not a divisor of NPROJS
    real,             parameter :: SMPD = 1.5, LAMBDA = 1.0e-3
    real,             parameter :: CROP_SMPD = SMPD * real(BOX) / real(CROP_BOX)
    real,             parameter :: MASS_LAMBDA_REL = 1.0e-3
    real,             parameter :: ADJOINT_RELTOL = 1.0e-5, NORMAL_OP_RELTOL = 1.0e-4
    real,             parameter :: RECON_CORR_THRES = 0.9, MONOTONIC_SLACK = 1.0e-3
    ! kernel interior tolerance: single-precision roundoff loosened for the kernel's
    ! shift-invariance approximation, not tuned on a reconstruction
    real,             parameter :: EPS_INTERIOR = 5.0e-2
    ! measure_kernel_scale returns 1.0 when calibrate_kernel's analytic padsc**2 is right
    real,             parameter :: KSCALE_TOL = 2.0e-2
    ! Fused and monolithic accumulation differ at single-precision roundoff;
    ! gate B and Khat directly and strictly. The preconditioner is a guarded reciprocal of D, and 20 Krylov
    ! recurrences amplify that harmless perturbation; give the final volume a
    ! separate 5e-4 bound so solver conditioning cannot hide an accumulator
    ! defect or manufacture a false failure.
    real,             parameter :: STREAM_ACCUM_RELTOL = 1.0e-6
    real,             parameter :: STREAM_SOLVE_RELTOL = 5.0e-4
    real,             parameter :: MASS_SCALE_RELTOL   = 5.0e-6
    real,             parameter :: MASS_SOLVE_RELTOL   = 5.0e-4
    real,             parameter :: CROP_RAW_RELTOL     = 2.0e-5
    integer,          parameter :: STREAM_ITS          = 20
    integer,          parameter :: KERNEL_COMPARE_ITS  = 8
    ! stage 13: the two backends' cropped observations differ only by the
    ! single-precision FFT round trips of the fused gridding route
    real,             parameter :: OBS_PARITY_RELTOL   = 1.0e-5
    ! stage 14: window band = mask3D_soft's soft ramp (values in
    ! [SUPPORT_BAND_LO,SUPPORT_BAND_HI]); repeated output-space warm starts
    ! must not apply the soft window more than once and shrink this band.
    ! SUPPORT_ITS stays at a production-sized budget: much longer solves of the small hard-masked
    ! operator reach a near-null mode of P(H+lambda I)P (PCG_STOP_INDEFINITE).
    real,             parameter :: SUPPORT_MSKRAD = real(BOX)/2.0 - 2.0
    integer,          parameter :: SUPPORT_ITS = 5, SUPPORT_NREP = 6
    real,             parameter :: SUPPORT_RTOL = 1.0e-4
    real,             parameter :: SUPPORT_BAND_LO = 0.15, SUPPORT_BAND_HI = 0.85
    real,             parameter :: SUPPORT_STABILITY_FRAC = 0.9
    real,             parameter :: CTRS(3,NBLOBS) = reshape([&
        &-5.0,-3.0, 2.0,&
        &4.0, 5.0,-3.0,&
        &0.0,-6.0,-5.0,&
        &3.0,-2.0, 6.0], [3,NBLOBS])
    real,             parameter :: SIGMAS(NBLOBS)    = [2.0, 2.5, 1.8, 2.2]
    real,             parameter :: AMPS(NBLOBS)      = [1.0, 0.8, 0.6, 0.5]
    real,             parameter :: KV = 300., CS = 2.7, FRACA = 0.1
    real,             parameter :: DFX_VALS(NCTF)    = [1.0, 1.5, 2.0, 2.5, 3.0]
    real,             parameter :: ASTIG_VALS(NCTF)  = [0.10, 0.15, 0.20, 0.12, 0.18]
    real,             parameter :: ANGAST_VALS(NCTF) = [0., 20., 40., 60., 80.]
    type(reconstructor_pcg) :: pcgop, pcg_reduce, pcg_crop, pcg_ml, pcg_c
    type(image)             :: ptcl_native, ptcl_work, obs_g, obs_p, pad_g, mskimg, wimg
    logical, allocatable    :: lmsk_native(:,:,:), lmsk_crop(:,:,:)
    real,    allocatable    :: rm_pad(:,:,:), obsg(:,:), obsp(:,:), proj2d(:,:,:), window(:,:,:)
    real,    allocatable    :: x_c(:,:,:)
    real    :: obs_scale, obs_err, rms_c, rms_first, rms_last, support_leak
    integer :: npad, irep
    type(oris)              :: projdirs, projdirs_exp, projdirs_crop
    type(ori)               :: e, e_exp
    type(ctfparams)         :: ctfparms
    type(sym)               :: c1sym, c2sym
    type(string)            :: raw_part1, raw_part2, raw_ml
    type(pcg_solver_outcome) :: solver_outcome
    real,    allocatable    :: phantom(:,:,:), p_probe(:,:,:), q_probe(:,:,:)
    real,    allocatable    :: hp(:,:,:), hq(:,:,:), hm(:,:,:), hk(:,:,:)
    real,    allocatable    :: recon(:,:,:), recon_str(:,:,:), rel_res_hist(:)
    real,    allocatable    :: recon_mf(:,:,:), recon_kernel(:,:,:)
    real,    allocatable    :: recon_mass(:,:,:), recon_mass_dup(:,:,:)
    real,    allocatable    :: p_crop(:,:,:), hm_crop(:,:,:), hk_crop(:,:,:)
    real,    allocatable    :: hm_ml(:,:,:), hk_ml(:,:,:), ml_prior_diag(:,:,:)
    real,    allocatable    :: xdiv(:,:,:), recon_on(:,:,:), recon_off(:,:,:)
    real,    allocatable    :: env(:,:,:), invenv(:,:,:)
    real,    allocatable    :: khat_a(:,:,:), khat_b(:,:,:), b_mono(:,:,:), b_str(:,:,:)
    real,    allocatable    :: qplane_re(:,:), qplane_im(:,:), sig2arr(:), sig2_2d(:,:)
    real,    allocatable    :: sig2_crop(:,:), draw_full(:,:,:), draw_crop(:,:,:), fsc_prior(:)
    complex, allocatable    :: gx_plane(:,:), qplane(:,:), mplane(:,:), wplane(:,:), Ti(:,:)
    complex, allocatable    :: adj_out(:,:,:), y_planes(:,:,:), y_exp(:,:,:), y_crop(:,:,:)
    complex, allocatable    :: braw_full(:,:,:), braw_crop(:,:,:)
    integer :: lims2(2,2), lims3(3,2), i, j, k, b, g, c, niters
    integer :: R, margin, lo, hi, ifrom, nb, nsym, nraw, nraw_total, ml_prior_npositive
    integer :: lims2_crop(2,2), lims3_crop(3,2), raw_lim, hraw, kraw, mraw
    real    :: ctr, dx, dy, dz, adjoint_err, corr, shift(2)
    real    :: err_all, err_int, err_max, den_all, den_int, kdiff, stream_err, rhs_err, kscale
    real    :: solution_err, solution_norm_ratio, energy_ratio
    real    :: corr_on, corr_off, env_ctr, env_edge, recip_err
    real    :: data_scale, data_scale_dup, lambda_eff, lambda_eff_dup
    real    :: mass_scale_err, mass_lambda_err, mass_solution_err
    real    :: crop_b_err, crop_d_err, crop_kernel_err, crop_factor
    real    :: ml_kernel_err, ml_prior_min, ml_prior_max
    real    :: ml_prior_positive_min, ml_prior_positive_max, ml_prior_to_khat_l1, ml_prior_to_khat_rms
    real(dp):: lhs, rhs, dp_p_hq, dp_hp_q, dp_p_hp
    real(dp):: crop_b_num, crop_b_den, crop_d_num, crop_d_den, ml_prior_energy
    logical :: all_ok
    all_ok = .true.

    ! ---- deterministic RNG seed: every probe below must be reproducible or a
    !      tolerance failure cannot be told from a different random draw ----
    call set_fixed_seed(42)

    ! ---- deterministic, asymmetric phantom (sum of off-centre Gaussian blobs).
    !      Asymmetric on purpose: a symmetric phantom hides orientation bugs. ----
    write(logfhandle,'(a)') '>>> TEST_PCG_RECON: building deterministic phantom'
    allocate(phantom(BOX,BOX,BOX), source=0.0)
    ctr = real(BOX)/2.0 + 0.5
    do k = 1,BOX
        do j = 1,BOX
            do i = 1,BOX
                do b = 1,NBLOBS
                    dx = real(i)-ctr-CTRS(1,b); dy = real(j)-ctr-CTRS(2,b); dz = real(k)-ctr-CTRS(3,b)
                    phantom(i,j,k) = phantom(i,j,k) + AMPS(b)*exp(-(dx*dx+dy*dy+dz*dz)/(2.0*SIGMAS(b)**2))
                end do
            end do
        end do
    end do

    call pcgop%new(BOX, SMPD, LAMBDA)
    ! Deapodization OFF for stages 1-7 -- see the inverse-crime note in the
    ! header. Stage 8 turns it back on and supplies envelope-free data.
    call pcgop%set_deapod(.false.)
    lims2 = pcgop%get_lims2()
    lims3 = pcgop%get_lims3()
    R     = lims2(1,2)
    allocate(sig2arr(0:R))
    do i = 0, R
        sig2arr(i) = 1.0 + 0.15*real(i)
    end do
    call e%new(.false.)

    ! ================= STAGE 1: adjoint identity, T = 1 =================
    write(logfhandle,'(a)') '>>> STAGE 1: adjoint dot-product identity (T = 1)'
    call e%set_euler([30.,55.,70.])
    allocate(gx_plane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    call pcgop%set_volume(phantom)
    call pcgop%forward_plane(e, gx_plane)
    allocate(qplane_re(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    allocate(qplane_im(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    call random_number(qplane_re); call random_number(qplane_im)
    allocate(qplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
    qplane = cmplx(qplane_re-0.5, qplane_im-0.5)
    lhs = sum(real(conjg(gx_plane)*qplane, dp))
    allocate(adj_out(lims3(1,1):lims3(1,2), lims3(2,1):lims3(2,2), lims3(3,1):lims3(3,2)), source=cmplx(0.,0.))
    call pcgop%adjoint_plane_add(qplane, e, adj_out)
    ! <x, G^dagger q> over the operator's own oversampled lattice; fourier_dot
    ! keeps that lattice and its Friedel packing an implementation detail
    rhs = pcgop%fourier_dot(adj_out)
    adjoint_err = real(abs(lhs-rhs) / max(1.0_dp, abs(lhs), abs(rhs)))
    write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    <Gx,q>=', real(lhs), ' <x,G^Tq>=', real(rhs), &
        &' rel_err=', adjoint_err
    if( adjoint_err > ADJOINT_RELTOL )then
        write(logfhandle,'(a)') '    FAIL: adjoint dot-product identity violated'
        all_ok = .false.
    else
        write(logfhandle,'(a)') '    PASS: adjoint dot-product identity holds'
    endif

    ! ============ STAGE 2: adjoint identity, weighted transfer ============
    ! The T = 1 case above is trivially self-adjoint and cannot exercise
    ! build_transfer; this one carries a real astigmatic CTF, a nonzero shift
    ! and a sigma2 profile.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 2: adjoint identity with T = C*S/sqrt(sigma2)'
        ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
        ctfparms%dfx = DFX_VALS(1); ctfparms%dfy = DFX_VALS(1)+ASTIG_VALS(1)
        ctfparms%angast = ANGAST_VALS(1); ctfparms%phshift = 0.
        shift = [1.7, -2.3]
        allocate(Ti(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        Ti = pcgop%build_transfer(ctfparms, shift, sig2arr)
        allocate(mplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        mplane = Ti * gx_plane
        lhs = sum(real(conjg(mplane)*qplane, dp))
        allocate(wplane(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2)))
        wplane  = conjg(Ti) * qplane
        adj_out = cmplx(0.,0.)
        call pcgop%adjoint_plane_add(wplane, e, adj_out)
        rhs = pcgop%fourier_dot(adj_out)
        adjoint_err = real(abs(lhs-rhs) / max(1.0_dp, abs(lhs), abs(rhs)))
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    <T*Gx,q>=', real(lhs), &
            &' <x,G^T(conjg(T)*q)>=', real(rhs), ' rel_err=', adjoint_err
        if( adjoint_err > ADJOINT_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: weighted adjoint identity violated (build_transfer)'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: weighted adjoint identity holds'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 2 SKIPPED: stage 1 failed'
    endif

    ! ---- heterogeneous selection used by stages 3-7 ----
    call projdirs%new(NPROJS, .false.)
    call projdirs%spiral
    do i = 1, NPROJS
        call projdirs%get_ori(i, e)
        g = mod(i-1, NCTF) + 1
        ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
        ctfparms%dfx = DFX_VALS(g); ctfparms%dfy = DFX_VALS(g)+ASTIG_VALS(g)
        ctfparms%angast = ANGAST_VALS(g); ctfparms%phshift = 0.
        call e%set_ctfvars(ctfparms)
        call e%set_shift([1.0+0.3*real(mod(i-1,5)), -1.0+0.2*real(mod(i-1,7))])
        call projdirs%set_ori(i, e)
    end do
    allocate(sig2_2d(0:R,NPROJS))
    do i = 1, NPROJS
        sig2_2d(:,i) = sig2arr
    end do

    ! ========== STAGE 3: normal operator, symmetry and PSD ==========
    ! CG requires both. Symmetry is checked by the dot-product identity across
    ! two independent random probes; PSD by a single quadratic form.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 3: normal-operator symmetry and positive-definiteness'
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        allocate(p_probe(BOX,BOX,BOX), q_probe(BOX,BOX,BOX))
        call random_number(p_probe); call random_number(q_probe)
        p_probe = p_probe - 0.5; q_probe = q_probe - 0.5
        hp = pcgop%apply_normal(p_probe)
        hq = pcgop%apply_normal(q_probe)
        dp_p_hq = pcgop%dot_real_volume(p_probe, hq)
        dp_hp_q = pcgop%dot_real_volume(hp, q_probe)
        dp_p_hp = pcgop%dot_real_volume(p_probe, hp)
        adjoint_err = real(abs(dp_p_hq-dp_hp_q) / max(1.0_dp, abs(dp_p_hq), abs(dp_hp_q)))
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    dot(p,Hq)=', real(dp_p_hq), &
            &' dot(Hp,q)=', real(dp_hp_q), ' rel_err=', adjoint_err
        if( adjoint_err > NORMAL_OP_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: normal operator is not symmetric'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: normal operator is symmetric'
        endif
        write(logfhandle,'(a,es14.6)') '    dot(p,Hp)=', real(dp_p_hp)
        if( dp_p_hp <= 0.0_dp )then
            write(logfhandle,'(a)') '    FAIL: normal operator is not positive-definite'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: normal operator is positive-definite'
        endif
        ! P H P with the support mask ACTIVE -- the operator production solves
        ! with (see set_window_sphere: x = P u gives (P H P) u = P b). The soft edge
        ! makes P non-idempotent, so a one-sided or asymmetric mask application
        ! breaks the dot-product identity where the unmasked check cannot see it.
        write(logfhandle,'(a)') '>>> STAGE 3b: masked-operator (P H P) symmetry and positive-definiteness'
        call pcgop%set_window_sphere(real(BOX)/3.0)
        hp = pcgop%apply_normal(p_probe)
        hq = pcgop%apply_normal(q_probe)
        dp_p_hq = pcgop%dot_real_volume(p_probe, hq)
        dp_hp_q = pcgop%dot_real_volume(hp, q_probe)
        dp_p_hp = pcgop%dot_real_volume(p_probe, hp)
        adjoint_err = real(abs(dp_p_hq-dp_hp_q) / max(1.0_dp, abs(dp_p_hq), abs(dp_hp_q)))
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    dot(p,PHPq)=', real(dp_p_hq), &
            &' dot(PHPp,q)=', real(dp_hp_q), ' rel_err=', adjoint_err
        if( adjoint_err > NORMAL_OP_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: masked normal operator is not symmetric'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: masked normal operator is symmetric'
        endif
        ! strict positivity holds because H carries the lambda ridge:
        ! p.(PHP)p = (Pp).H(Pp) >= lambda*|Pp|^2 > 0 for a random probe
        write(logfhandle,'(a,es14.6)') '    dot(p,PHPp)=', real(dp_p_hp)
        if( dp_p_hp <= 0.0_dp )then
            write(logfhandle,'(a)') '    FAIL: masked normal operator is not positive-definite'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: masked normal operator is positive-definite'
        endif
        call pcgop%set_window_sphere(0.0)   ! stages 4+ assert on the unmasked operator
    else
        write(logfhandle,'(a)') '>>> STAGE 3 SKIPPED: an earlier stage failed'
    endif

    ! ============== STAGE 4: heterogeneous synthetic recovery ==============
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 4: synthetic recovery through CTF + shift + sigma'
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS))
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        end do
        ! build_operators here rather than a bare solve, so stage 7 can compare
        ! the two accumulation paths with identical preconditioner state
        call pcgop%build_operators(.false.)
        allocate(recon(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon, maxits=40, rtol=1.0e-3, &
            &rel_res_hist=rel_res_hist, niters=niters)
        write(logfhandle,'(a,i0,a)') '    PCG ran ', niters, ' iterations'
        do i = 2, niters
            if( rel_res_hist(i) > rel_res_hist(i-1) + MONOTONIC_SLACK )then
                write(logfhandle,'(a,i0)') '    FAIL: relative residual increased at iteration ', i
                all_ok = .false.
            endif
        end do
        corr = corr_of(recon, phantom)
        write(logfhandle,'(a,f8.5)') '    recon-vs-phantom correlation = ', corr
        if( corr < RECON_CORR_THRES )then
            write(logfhandle,'(a)') '    FAIL: recovered volume does not correlate with the known phantom'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: recovered volume correlates with the known phantom'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 4 SKIPPED: an earlier stage failed'
    endif

    ! ===== STAGE 5: kernel vs matrix-free operator and solve baseline =====
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 5: kernelized vs matrix-free normal operator'
        call pcgop%build_kernel
        hm = pcgop%apply_normal_matrixfree(p_probe)
        hk = pcgop%apply_normal_kernel(p_probe)
        den_all = max(1.0, sqrt(sum(hm*hm)))
        err_all = sqrt(sum((hk-hm)**2)) / den_all
        ! interior margin at least the KB support, so the shift-invariance
        ! approximation is not judged on the boundary where it is known to differ
        margin  = 4
        lo      = margin + 1
        hi      = BOX - margin
        den_int = max(1.0, sqrt(sum(hm(lo:hi,lo:hi,lo:hi)**2)))
        err_int = sqrt(sum((hk(lo:hi,lo:hi,lo:hi)-hm(lo:hi,lo:hi,lo:hi))**2)) / den_int
        err_max = maxval(abs(hk-hm)) / max(1.0, maxval(abs(hm)))
        write(logfhandle,'(a,es14.6)') '    relative error, all voxels = ', err_all
        write(logfhandle,'(a,es14.6)') '    relative error, interior   = ', err_int
        write(logfhandle,'(a,es14.6)') '    relative maximum error     = ', err_max
        if( err_int > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: kernelized operator disagrees with the matrix-free reference'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernelized operator agrees with the reference in the interior'
        endif
        ! A scale error of a few percent would pass the interior check above but
        ! is exactly what a drift in the analytic constant would look like, so
        ! assert on it directly.
        kscale = pcgop%measure_kernel_scale()
        write(logfhandle,'(a,f10.6)') '    measured kernel scale (1.0 = analytic padsc**2 exact) = ', kscale
        if( abs(kscale - 1.0) > KSCALE_TOL )then
            write(logfhandle,'(a)') '    FAIL: kernel scale has drifted from the analytic constant'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: analytic kernel scale is correct'
        endif
        dp_p_hq = pcgop%dot_real_volume(p_probe, hm)
        dp_hp_q = pcgop%dot_real_volume(p_probe, hk)
        energy_ratio = real(dp_hp_q / max(abs(dp_p_hq), epsilon(1.0_dp)))
        write(logfhandle,'(a,es14.6)') '    kernel/matrix-free probe energy ratio = ', energy_ratio

        ! Fixed-iteration solve comparison: reported, not gated.
        ! Same RHS, initial state and iteration count prevent convergence-stop
        ! jitter from masquerading as an operator difference.
        allocate(recon_mf(BOX,BOX,BOX), recon_kernel(BOX,BOX,BOX), source=0.0)
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        call pcgop%solve(y_planes, recon_mf, maxits=KERNEL_COMPARE_ITS, rtol=0.0)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call pcgop%solve(y_planes, recon_kernel, maxits=KERNEL_COMPARE_ITS, rtol=0.0)
        solution_err = sqrt(sum((recon_kernel-recon_mf)**2)) / &
            &max(1.0, sqrt(sum(recon_mf*recon_mf)))
        solution_norm_ratio = sqrt(sum(recon_kernel*recon_kernel)) / &
            &max(1.0, sqrt(sum(recon_mf*recon_mf)))
        write(logfhandle,'(a,i0,a,es14.6)') '    fixed-', KERNEL_COMPARE_ITS, &
            &'-iteration solution rel_err = ', solution_err
        write(logfhandle,'(a,es14.6)') '    fixed-iteration solution norm ratio = ', solution_norm_ratio
        deallocate(recon_mf, recon_kernel)
    else
        write(logfhandle,'(a)') '>>> STAGE 5 SKIPPED: an earlier stage failed'
    endif

    ! ========= STAGE 6: kernel invariants and preconditioner =========
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 6: kernel invariants and preconditioner'
        khat_a = pcgop%apply_normal_kernel(p_probe)
        ! changing ONLY the shifts must leave the normal operator unchanged: the
        ! shift is a unit-modulus phase and cancels in |T|^2
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call e%set_shift([-2.0+0.11*real(i), 3.0-0.07*real(i)])
            call projdirs%set_ori(i, e)
        end do
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_b = pcgop%apply_normal_kernel(p_probe)
        kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,es14.6)') '    relative change after shift-only edit = ', kdiff
        if( kdiff > 1.0e-5 )then
            write(logfhandle,'(a)') '    FAIL: kernel changed when only shifts changed'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernel is shift-invariant'
        endif
        ! changing a CTF MUST change it
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            ctfparms = e%get_ctfvars()
            ctfparms%dfx = ctfparms%dfx + 0.7
            ctfparms%dfy = ctfparms%dfy + 0.7
            call e%set_ctfvars(ctfparms)
            call projdirs%set_ori(i, e)
        end do
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_b = pcgop%apply_normal_kernel(p_probe)
        kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,es14.6)') '    relative change after CTF edit        = ', kdiff
        if( kdiff < 1.0e-4 )then
            write(logfhandle,'(a)') '    FAIL: kernel did not change when the CTF changed'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: kernel tracks the CTF'
        endif
        call pcgop%set_op_mode(PCG_OP_MATRIXFREE)
        call pcgop%build_precond
        hm = pcgop%apply_normal_matrixfree(p_probe)
        if( any(hm /= hm) )then
            write(logfhandle,'(a)') '    FAIL: operator produced non-finite values after build_precond'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: preconditioner built, operator still finite'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 6 SKIPPED: an earlier stage failed'
    endif

    ! ================== STAGE 7: streaming accumulation ==================
    ! begin_accum/accumulate_batch/end_accum/solve_accum must reproduce what
    ! solve() produces from all planes at once. The accumulated B and D-derived
    ! kernel must agree near single-precision roundoff; the fixed-step solution
    ! has its own bound because the reciprocal preconditioner and CG recurrence
    ! amplify that input perturbation.
    ! BATCHSZ is deliberately not a divisor of NPROJS so the short final batch
    ! is exercised.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 7: streaming accumulation vs monolithic solve'
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        ! regenerate observations: stage 6 edited the shifts and CTFs
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        end do
        ! BOTH paths run a FIXED number of iterations with the tolerance
        ! disabled. Comparing two tolerance-stopped solves is meaningless:
        ! round-off in the accumulators can put the stopping decision one
        ! iteration apart, and a single CG step near convergence moves the
        ! solution by O(rtol) -- orders of magnitude above anything worth
        ! asserting. Forcing identical steps leaves accumulator round-off as the
        ! only difference between the two paths, which is what is being tested.
        ! reference path: one shot
        call pcgop%build_operators(.false.)
        recon = 0.0
        call pcgop%solve(y_planes, recon, maxits=STREAM_ITS, rtol=0.0, &
            &niters=niters, outcome=solver_outcome)
        call pcgop%get_rhs(b_mono)
        write(logfhandle,'(a,i0,a)') '    monolithic: ', niters, ' iterations'
        if( trim(solver_outcome%stop_reason) /= 'fixed_iterations' .or. &
            &solver_outcome%iteration_count /= STREAM_ITS .or. solver_outcome%converged )then
            write(logfhandle,'(a)') '    FAIL: fixed-iteration solver outcome is inconsistent'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: fixed-iteration solver outcome is explicit'
        endif
        ! Two streaming runs, so a failure localizes itself rather than just
        ! reporting a number: ONE batch covering everything isolates the
        ! begin/end_accum machinery, and many batches adds the per-batch
        ! indexing on top. If single-batch agrees and multi-batch does not, the
        ! fault is in the batch range handling; if neither agrees, it is in the
        ! accumulate/fold path shared by both.
        allocate(recon_str(BOX,BOX,BOX))
        do j = 1, 2
            if( j == 1 )then
                nb = NPROJS                     ! single batch
            else
                nb = BATCHSZ                    ! many batches, short final one
            endif
            call pcgop%begin_accum
            do ifrom = 1, NPROJS, nb
                k = min(nb, NPROJS - ifrom + 1)
                call pcgop%accumulate_batch(y_planes(:,:,ifrom:ifrom+k-1), k, ifrom)
            end do
            call pcgop%end_accum(.false.)
            ! compare the RHS before the solve: if b already differs the fault is
            ! in fused B accumulation or the fold; if b agrees but the solutions
            ! do not, it is in the fused D accumulator feeding the preconditioner
            call pcgop%get_rhs(b_str)
            rhs_err = sqrt(sum((b_str-b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
            recon_str = 0.0
            call pcgop%solve_accum(recon_str, maxits=STREAM_ITS, rtol=0.0, niters=niters)
            stream_err = sqrt(sum((recon_str-recon)**2)) / max(1.0, sqrt(sum(recon*recon)))
            write(logfhandle,'(a,i3,a,es14.6,a,es14.6)') '    batch size ', nb, &
                &': rel_err(b) = ', rhs_err, '   rel_err(x) = ', stream_err
            if( rhs_err > STREAM_ACCUM_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: fused RHS does not reproduce monolithic accumulation'
                all_ok = .false.
            else if( stream_err > STREAM_SOLVE_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: fused accumulation changes the fixed-step solve beyond tolerance'
                all_ok = .false.
            else
                write(logfhandle,'(a)') '    PASS: streaming accumulation reproduces the monolithic solve'
            endif
        end do
        ! Production always uses the kernel operator, which reaches it through
        ! end_accum(.true.) -- a branch the runs above never touch, since they
        ! pass .false. Deriving Khat from a streamed accumulator must match
        ! deriving it from a monolithic one, or the default path is untested.
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%build_kernel
        khat_a = pcgop%apply_normal_kernel(p_probe)
        do j = 1, 2
            if( j == 1 )then
                nb = NPROJS
            else
                nb = BATCHSZ
            endif
            call pcgop%begin_accum
            do ifrom = 1, NPROJS, nb
                k = min(nb, NPROJS - ifrom + 1)
                call pcgop%accumulate_batch(y_planes(:,:,ifrom:ifrom+k-1), k, ifrom)
            end do
            call pcgop%end_accum(.true.)
            khat_b = pcgop%apply_normal_kernel(p_probe)
            kdiff  = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
            write(logfhandle,'(a,i3,a,es14.6)') '    kernel batch size ', nb, &
                &': streamed vs monolithic rel_err = ', kdiff
            if( kdiff > STREAM_ACCUM_RELTOL )then
                write(logfhandle,'(a)') '    FAIL: kernel built from a streamed accumulator differs'
                all_ok = .false.
            else
                write(logfhandle,'(a)') '    PASS: kernel is identical either way'
            endif
        end do
        ! Distributed contract: workers publish raw, unfolded B and D. The
        ! master adds artifacts in ascending part order and only then folds and
        ! finalizes. Split the same particle sequence into two artifacts so the
        ! gate covers serialization, fixed-order reduction and master-only
        ! finalization without involving a scheduler.
        raw_part1 = 'test_pcg_raw_part1.dat'
        raw_part2 = 'test_pcg_raw_part2.dat'
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes(:,:,1:NPROJS/2), NPROJS/2, 1)
        call pcgop%write_raw_accum(raw_part1, 1, 0, 1, 2, NPROJS/2, 'pcg_recon_test_v1')
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes(:,:,NPROJS/2+1:NPROJS), NPROJS/2, NPROJS/2+1)
        call pcgop%write_raw_accum(raw_part2, 1, 0, 2, 2, NPROJS/2, 'pcg_recon_test_v1')
        call pcgop%end_accum(.false.)
        call pcg_reduce%new(BOX, SMPD, LAMBDA)
        call pcg_reduce%set_deapod(.false.)
        call pcg_reduce%begin_reduction
        nraw_total = 0
        call pcg_reduce%add_raw_accum(raw_part1, 1, 0, 1, 2, 'pcg_recon_test_v1', nraw)
        nraw_total = nraw_total + nraw
        call pcg_reduce%add_raw_accum(raw_part2, 1, 0, 2, 2, 'pcg_recon_test_v1', nraw)
        nraw_total = nraw_total + nraw
        call pcg_reduce%end_accum(.true.)
        call pcg_reduce%get_rhs(b_str)
        khat_b = pcg_reduce%apply_normal_kernel(p_probe)
        rhs_err = sqrt(sum((b_str-b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
        kdiff   = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        write(logfhandle,'(a,i0,a,es14.6,a,es14.6)') '    raw reduction particles ', nraw_total, &
            &': rel_err(b) = ', rhs_err, '   rel_err(kernel) = ', kdiff
        if( nraw_total /= NPROJS .or. rhs_err > STREAM_ACCUM_RELTOL .or. kdiff > STREAM_ACCUM_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: raw fixed-order reduction differs from monolithic accumulation'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: raw fixed-order reduction reproduces monolithic accumulation'
        endif
        call pcg_reduce%kill
        call del_file(raw_part1)
        call del_file(raw_part2)
        call raw_part1%kill
        call raw_part2%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 7 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 8: deapodization, the honest stage ============
    ! It flips deapod state and rebuilds the selection.
    ! forward_plane(x/env) = FT[env . (x/env)] = FT[x], i.e. the true
    ! envelope-free central section -- no inverse crime.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 8: deapodization against envelope-free observations'
        env    = pcgop%get_env()
        invenv = pcgop%get_invenv()
        c         = BOX/2 + 1
        env_ctr   = env(c,c,c)
        env_edge  = env(1,c,c)
        recip_err = maxval(abs(env*invenv - 1.0))
        write(logfhandle,'(a,f10.6)')  '    envelope at centre       = ', env_ctr
        write(logfhandle,'(a,f10.6)')  '    envelope at edge (1,c,c) = ', env_edge
        write(logfhandle,'(a,es14.6)') '    max|env*invenv - 1|      = ', recip_err
        if( abs(env_ctr - 1.0) > 1.0e-5 )then
            write(logfhandle,'(a)') '    FAIL: envelope is not unity at the box centre'
            all_ok = .false.
        else if( env_edge >= env_ctr )then
            write(logfhandle,'(a)') '    FAIL: envelope does not decay away from the centre'
            all_ok = .false.
        else if( recip_err > 1.0e-4 )then
            write(logfhandle,'(a)') '    FAIL: invenv is not the reciprocal of env'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: envelope measured, normalized and invertible'
        endif
    endif
    if( all_ok )then
        allocate(xdiv(BOX,BOX,BOX))
        xdiv = phantom * invenv
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call pcgop%set_volume(xdiv)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        call pcgop%set_deapod(.true.)
        call pcgop%prep_particles(projdirs)
        allocate(recon_on(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon_on, maxits=40, rtol=1.0e-3, niters=niters)
        corr_on = corr_of(recon_on, phantom)
        write(logfhandle,'(a,i0,a,f9.5)') '    deapod ON : ', niters, ' iters, corr = ', corr_on
        call pcgop%set_deapod(.false.)
        call pcgop%prep_particles(projdirs)
        allocate(recon_off(BOX,BOX,BOX), source=0.0)
        call pcgop%solve(y_planes, recon_off, maxits=40, rtol=1.0e-3, niters=niters)
        corr_off = corr_of(recon_off, phantom)
        write(logfhandle,'(a,i0,a,f9.5)') '    deapod OFF: ', niters, ' iters, corr = ', corr_off
        if( corr_on <= corr_off )then
            write(logfhandle,'(a)') '    FAIL: deapodization did not improve recovery of the true volume'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: deapodization recovers the true volume better'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 8 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 9: symmetry by coordinate replication ============
    ! Replication over c2 must build the same Khat and b as the c2-expanded particle set at c1 (policy
    ! sec. 6). Compared at operator level, not through a solve: the kernelized operator is only nearly SPD.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 9: symmetry replication == particle-set expansion (c2, kernel)'
        call c1sym%new('c1')
        call c2sym%new('c2')
        nsym = c2sym%get_nsym()
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        call pcgop%set_deapod(.true.)   ! stage 8 leaves it OFF; test the shipped config
        ! fresh, well-conditioned data: spiral projections of the phantom, T = 1
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        if( allocated(y_planes) ) deallocate(y_planes)
        allocate(y_planes(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS))
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        ! Path A: in-operator replication over c2 on NPROJS particles
        call pcgop%set_sym(c2sym)
        call pcgop%prep_particles(projdirs)
        call pcgop%build_operators(.true.)
        khat_a = pcgop%apply_normal_kernel(p_probe)
        b_mono = pcgop%apply_adjoint_all(y_planes)
        ! Path B: explicit c1 build of the c2-expanded particle set. The mate
        ! orientation is composed through sym%apply -- production's own symmetry
        ! adaptor, the one commander_reconstruct3D reaches through insert_fplane
        ! -- NOT by repeating the reconstructor's internal matmul(R_i, S_g). That
        ! is the point: composing both paths the same way would make this stage
        ! agree under a transposed or reversed composition too, and prove
        ! nothing. Going through sym%apply makes the stage an assertion that the
        ! operator's convention IS production's.
        call projdirs_exp%new(NPROJS*nsym, .false.)
        allocate(y_exp(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), NPROJS*nsym))
        c = 0
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            do g = 1, nsym
                c = c + 1
                call c2sym%apply(e, g, e_exp)
                call projdirs_exp%set_ori(c, e_exp)
                y_exp(:,:,c) = y_planes(:,:,i)
            end do
        end do
        call pcgop%set_sym(c1sym)          ! back to no replication
        call pcgop%prep_particles(projdirs_exp)
        call pcgop%build_operators(.true.)
        khat_b = pcgop%apply_normal_kernel(p_probe)
        b_str  = pcgop%apply_adjoint_all(y_exp)
        kdiff      = sqrt(sum((khat_b-khat_a)**2)) / max(1.0, sqrt(sum(khat_a*khat_a)))
        stream_err = sqrt(sum((b_str -b_mono)**2)) / max(1.0, sqrt(sum(b_mono*b_mono)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    rel_err(Khat) = ', kdiff, '   rel_err(b) = ', stream_err
        if( kdiff > 5.0e-3 .or. stream_err > 5.0e-3 )then
            write(logfhandle,'(a)') '    FAIL: symmetry replication does not match the expanded-set build'
            all_ok = .false.
        else if( sqrt(sum(khat_a*khat_a)) <= TINY )then
            write(logfhandle,'(a)') '    FAIL: replicated kernel is trivially zero'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: symmetry replication equals the expanded-set build'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 9 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 10: lambda scaling with effective data mass ============
    ! If every weighted observation is duplicated, B and H_data both double.
    ! A fixed absolute lambda would then become half as strong. Relative lambda
    ! must double with D so the complete normal system is scaled uniformly and
    ! a fixed-iteration PCG solve remains invariant.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 10: lambda scaling with effective data mass'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call pcgop%set_lambda_relative(MASS_LAMBDA_REL)
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, y_planes(:,:,i))
        end do
        call pcgop%prep_particles(projdirs)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes, NPROJS, 1)
        call pcgop%end_accum(.true.)
        call pcgop%set_op_mode(PCG_OP_KERNEL)
        data_scale = pcgop%get_data_scale()
        lambda_eff = pcgop%get_effective_lambda()
        allocate(recon_mass(BOX,BOX,BOX), source=0.0)
        call pcgop%solve_accum(recon_mass, maxits=2, rtol=0.0)

        call projdirs_exp%kill
        call projdirs_exp%new(2*NPROJS, .false.)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call projdirs_exp%set_ori(i, e)
            call projdirs_exp%set_ori(NPROJS+i, e)
        end do
        if( allocated(y_exp) ) deallocate(y_exp)
        allocate(y_exp(lims2(1,1):lims2(1,2), lims2(2,1):lims2(2,2), 2*NPROJS))
        y_exp(:,:,1:NPROJS)          = y_planes
        y_exp(:,:,NPROJS+1:2*NPROJS) = y_planes
        ! Reach the doubled-data solve through the production distributed
        ! boundary: two raw worker artifacts, fixed-order master reduction,
        ! then relative-lambda derivation from the reduced D.
        raw_part1 = 'test_pcg_mass_part1.dat'
        raw_part2 = 'test_pcg_mass_part2.dat'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call pcgop%prep_particles(projdirs_exp)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_exp(:,:,1:NPROJS), NPROJS, 1)
        call pcgop%write_raw_accum(raw_part1, 1, 0, 1, 2, NPROJS, 'pcg_lambda_mass_v1')
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_exp(:,:,NPROJS+1:2*NPROJS), NPROJS, NPROJS+1)
        call pcgop%write_raw_accum(raw_part2, 1, 0, 2, 2, NPROJS, 'pcg_lambda_mass_v1')
        call pcgop%kill
        call pcg_reduce%new(BOX, SMPD, 0.0)
        call pcg_reduce%set_deapod(.false.)
        call pcg_reduce%set_lambda_relative(MASS_LAMBDA_REL)
        call pcg_reduce%begin_reduction
        call pcg_reduce%add_raw_accum(raw_part1, 1, 0, 1, 2, 'pcg_lambda_mass_v1', nraw)
        call pcg_reduce%add_raw_accum(raw_part2, 1, 0, 2, 2, 'pcg_lambda_mass_v1', nraw)
        call pcg_reduce%end_accum(.true.)
        call pcg_reduce%set_op_mode(PCG_OP_KERNEL)
        data_scale_dup = pcg_reduce%get_data_scale()
        lambda_eff_dup = pcg_reduce%get_effective_lambda()
        allocate(recon_mass_dup(BOX,BOX,BOX), source=0.0)
        call pcg_reduce%solve_accum(recon_mass_dup, maxits=2, rtol=0.0)

        mass_scale_err = abs(data_scale_dup - 2.0*data_scale) / max(TINY, 2.0*data_scale)
        mass_lambda_err = abs(lambda_eff_dup - 2.0*lambda_eff) / max(TINY, 2.0*lambda_eff)
        mass_solution_err = sqrt(sum((recon_mass_dup-recon_mass)**2)) / &
            &max(1.0, sqrt(sum(recon_mass*recon_mass)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    data scale: original = ', data_scale, &
            &' duplicated = ', data_scale_dup
        write(logfhandle,'(a,es14.6,a,es14.6)') '    lambda_eff: original = ', lambda_eff, &
            &' duplicated = ', lambda_eff_dup
        write(logfhandle,'(a,es14.6,a,es14.6,a,es14.6)') '    rel_err(scale x2) = ', mass_scale_err, &
            &' rel_err(lambda x2) = ', mass_lambda_err, ' rel_err(solution) = ', mass_solution_err
        if( mass_scale_err > MASS_SCALE_RELTOL .or. mass_lambda_err > MASS_SCALE_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: effective lambda does not scale with duplicated data'
            all_ok = .false.
        else if( mass_solution_err > MASS_SOLVE_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: relative lambda changes the solution under data duplication'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: relative lambda preserves the solution across data mass'
        endif
        call del_file(raw_part1)
        call del_file(raw_part2)
        call raw_part1%kill
        call raw_part2%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 10 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 11: cropped accumulator and kernel invariants ============
    ! A cropped reconstruction covers the same physical extent with fewer
    ! Fourier samples. Native and cropped operators must therefore deposit the
    ! same raw sufficient statistics inside the common band, provided the crop
    ! boundary and the KB stencil are excluded. The cropped Khat must also
    ! retain the matrix-free equivalence already required at the native box.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 11: full vs cropped raw accumulators and cropped kernel'
        call pcgop%kill
        call pcgop%new(BOX, SMPD, 0.0)
        call pcgop%set_deapod(.false.)
        call projdirs%kill
        call projdirs%new(NPROJS, .false.)
        call projdirs%spiral
        call projdirs_crop%new(NPROJS, .false.)
        crop_factor = real(CROP_BOX) / real(BOX)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            g = mod(i-1, NCTF) + 1
            ctfparms%smpd = SMPD; ctfparms%kv = KV; ctfparms%cs = CS; ctfparms%fraca = FRACA
            ctfparms%dfx = DFX_VALS(g); ctfparms%dfy = DFX_VALS(g)+ASTIG_VALS(g)
            ctfparms%angast = ANGAST_VALS(g); ctfparms%phshift = 0.
            shift = [1.0+0.3*real(mod(i-1,5)), -1.0+0.2*real(mod(i-1,7))]
            call e%set_ctfvars(ctfparms)
            call e%set_shift(shift)
            call projdirs%set_ori(i, e)
            ctfparms%smpd = CROP_SMPD
            call e%set_ctfvars(ctfparms)
            call e%set_shift(shift*crop_factor)
            call projdirs_crop%set_ori(i, e)
        enddo
        call pcgop%set_volume(phantom)
        do i = 1, NPROJS
            call projdirs%get_ori(i, e)
            call pcgop%forward_plane(e, gx_plane)
            Ti = pcgop%build_transfer(e%get_ctfvars(), e%get_2Dshift(), sig2arr)
            y_planes(:,:,i) = Ti * gx_plane
        enddo
        call pcgop%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcgop%begin_accum
        call pcgop%accumulate_batch(y_planes, NPROJS, 1)
        call pcgop%get_raw_accum(braw_full, draw_full)

        call pcg_crop%new(CROP_BOX, CROP_SMPD, 0.0)
        call pcg_crop%set_deapod(.true.)
        lims2_crop = pcg_crop%get_lims2()
        lims3_crop = pcg_crop%get_lims3()
        allocate(sig2_crop(0:CROP_BOX/2,NPROJS))
        do i = 1, NPROJS
            sig2_crop(:,i) = sig2arr(0:CROP_BOX/2)
        enddo
        allocate(y_crop(lims2_crop(1,1):lims2_crop(1,2), &
                       &lims2_crop(2,1):lims2_crop(2,2), NPROJS))
        y_crop = y_planes(lims2_crop(1,1):lims2_crop(1,2), &
            &lims2_crop(2,1):lims2_crop(2,2),:)
        call pcg_crop%prep_particles(projdirs_crop, use_ctf=.true., sig2=sig2_crop)
        call pcg_crop%begin_accum
        call pcg_crop%accumulate_batch(y_crop, NPROJS, 1)
        call pcg_crop%get_raw_accum(braw_crop, draw_crop)
        raw_ml = 'test_pcg_ml_prior.dat'
        call pcg_crop%write_raw_accum(raw_ml, 1, 0, 1, 1, NPROJS, 'pcg_ml_prior_v1')

        raw_lim    = CROP_BOX - 5
        crop_b_num = 0.0_dp
        crop_b_den = 0.0_dp
        crop_d_num = 0.0_dp
        crop_d_den = 0.0_dp
        do mraw = lims3_crop(3,1), lims3_crop(3,2)
            do kraw = lims3_crop(2,1), lims3_crop(2,2)
                do hraw = lims3_crop(1,1), lims3_crop(1,2)
                    if( hraw*hraw + kraw*kraw + mraw*mraw > raw_lim*raw_lim ) cycle
                    crop_b_num = crop_b_num + real(abs(braw_crop(hraw,kraw,mraw) - &
                        &braw_full(hraw,kraw,mraw))**2,dp)
                    crop_b_den = crop_b_den + real(abs(braw_full(hraw,kraw,mraw))**2,dp)
                    crop_d_num = crop_d_num + real((draw_crop(hraw,kraw,mraw) - &
                        &draw_full(hraw,kraw,mraw))**2,dp)
                    crop_d_den = crop_d_den + real(draw_full(hraw,kraw,mraw)**2,dp)
                enddo
            enddo
        enddo
        crop_b_err = real(sqrt(crop_b_num / max(1.0_dp,crop_b_den)))
        crop_d_err = real(sqrt(crop_d_num / max(1.0_dp,crop_d_den)))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    common-band rel_err(B) = ', crop_b_err, &
            &' rel_err(D) = ', crop_d_err
        if( crop_b_err > CROP_RAW_RELTOL .or. crop_d_err > CROP_RAW_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: cropped raw accumulators disagree with the full-box common band'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: cropped raw accumulators preserve the full-box common band'
        endif

        call pcg_crop%end_accum(.true.)
        allocate(p_crop(CROP_BOX,CROP_BOX,CROP_BOX))
        call random_number(p_crop)
        p_crop = p_crop - 0.5
        hm_crop = pcg_crop%apply_normal_matrixfree(p_crop)
        hk_crop = pcg_crop%apply_normal_kernel(p_crop)
        margin = 4
        lo = margin + 1
        hi = CROP_BOX - margin
        crop_kernel_err = sqrt(sum((hk_crop(lo:hi,lo:hi,lo:hi) - &
            &hm_crop(lo:hi,lo:hi,lo:hi))**2)) / &
            &max(1.0,sqrt(sum(hm_crop(lo:hi,lo:hi,lo:hi)**2)))
        write(logfhandle,'(a,es14.6)') '    cropped kernel interior rel_err = ', crop_kernel_err
        if( crop_kernel_err > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: cropped kernel disagrees with the matrix-free reference'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: cropped kernel agrees with the matrix-free reference'
        endif
    else
        write(logfhandle,'(a)') '>>> STAGE 11 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 12: FSC/SSNR ML prior operator ============
    ! Replay the cropped raw artifact so the prior is derived only after D is
    ! available, as the production two-map path does after the base
    ! independent-half FSC has been measured.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 12: FSC/SSNR ML prior positivity and operator parity'
        call pcg_ml%new(CROP_BOX, CROP_SMPD, 0.0)
        call pcg_ml%set_deapod(.true.)
        call pcg_ml%prep_particles(projdirs_crop, use_ctf=.true., sig2=sig2_crop)
        call pcg_ml%begin_reduction
        call pcg_ml%add_raw_accum(raw_ml, 1, 0, 1, 1, 'pcg_ml_prior_v1', nraw)
        allocate(fsc_prior(CROP_BOX/2), source=0.5)
        do i = 1, size(fsc_prior)
            fsc_prior(i) = max(0.05, 0.9 - 0.1*real(i-1))
        enddo
        call pcg_ml%set_ml_prior(fsc_prior, 1.0, 100.0)
        call pcg_ml%end_accum(.true.)
        call pcg_ml%get_ml_prior(ml_prior_diag)
        call pcg_ml%get_ml_prior_stats(ml_prior_npositive, ml_prior_positive_min, ml_prior_positive_max, &
            &ml_prior_to_khat_l1, ml_prior_to_khat_rms)
        ml_prior_min = minval(ml_prior_diag)
        ml_prior_max = maxval(ml_prior_diag)
        write(logfhandle,'(a,es14.6,a,es14.6)') '    prior diagonal min = ', ml_prior_min, &
            &' max = ', ml_prior_max
        if( ml_prior_min < 0.0 .or. ml_prior_max <= 0.0 )then
            write(logfhandle,'(a)') '    FAIL: ML prior is not a nonzero positive-semidefinite diagonal'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior is a nonzero positive-semidefinite diagonal'
        endif
        write(logfhandle,'(a,i0,4(a,es14.6))') '    positive bins = ', ml_prior_npositive, &
            &' min+ = ', ml_prior_positive_min, ' max+ = ', ml_prior_positive_max, &
            &' prior/Khat L1 = ', ml_prior_to_khat_l1, ' RMS = ', ml_prior_to_khat_rms
        if( ml_prior_npositive /= count(ml_prior_diag > 0.0) .or. &
            &abs(ml_prior_positive_min-minval(ml_prior_diag, mask=ml_prior_diag > 0.0)) > 0.0 .or. &
            &abs(ml_prior_positive_max-ml_prior_max) > 0.0 .or. &
            &ml_prior_to_khat_l1 <= 0.0 .or. ml_prior_to_khat_rms <= 0.0 )then
            write(logfhandle,'(a)') '    FAIL: ML prior summary diagnostics are inconsistent'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior summary diagnostics are consistent'
        endif
        hm_ml = pcg_ml%apply_normal_matrixfree(p_crop)
        hk_ml = pcg_ml%apply_normal_kernel(p_crop)
        margin = 4
        lo = margin + 1
        hi = CROP_BOX - margin
        ml_kernel_err = sqrt(sum((hk_ml(lo:hi,lo:hi,lo:hi) - &
            &hm_ml(lo:hi,lo:hi,lo:hi))**2)) / &
            &max(1.0,sqrt(sum(hm_ml(lo:hi,lo:hi,lo:hi)**2)))
        ml_prior_energy = sum(real(p_crop,dp) * real(hk_ml-hk_crop,dp))
        write(logfhandle,'(a,es14.6,a,es14.6)') '    ML kernel interior rel_err = ', ml_kernel_err, &
            &' dot(p,P_tau*p) = ', real(ml_prior_energy)
        if( ml_kernel_err > EPS_INTERIOR )then
            write(logfhandle,'(a)') '    FAIL: ML kernel disagrees with the matrix-free reference'
            all_ok = .false.
        else if( ml_prior_energy <= 0.0_dp )then
            write(logfhandle,'(a)') '    FAIL: ML prior does not add positive quadratic energy'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: ML prior preserves operator parity and adds positive energy'
        endif
        call del_file(raw_ml)
        call raw_ml%kill
    else
        write(logfhandle,'(a)') '>>> STAGE 12 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 13: prepared-observation parity (gridding vs PCG) ============
    ! Cropping and edge tapering do not commute, so both backends prepare a
    ! cropped particle through prep_rec_observation. The gridding route is the
    ! fused crop/taper/pad/FFT chain of prep_imgs4rec; that routine leaves its
    ! cropped input in the tapered real-space state. Compare that input with
    ! the PCG preparation before backend-specific transforms or accumulation.
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 13: prepared-observation parity for box_crop < box'
        call ptcl_native%new([BOX,BOX,1], SMPD)
        call ptcl_work%new([BOX,BOX,1], SMPD)
        call obs_g%new([CROP_BOX,CROP_BOX,1], CROP_SMPD)
        call obs_p%new([CROP_BOX,CROP_BOX,1], CROP_SMPD)
        npad = OSMPL_PAD_FAC * CROP_BOX
        call pad_g%new([npad,npad,1], CROP_SMPD)
        call mskimg%disc([BOX,BOX,1], SMPD, real(BOX)/3.0, lmsk_native)
        call mskimg%disc([CROP_BOX,CROP_BOX,1], CROP_SMPD, real(CROP_BOX)/3.0, lmsk_crop)
        call mskimg%kill
        ! deterministic particle: the phantom projected along z on a background ramp plus noise
        call ptcl_native%gauran(0., 0.3)
        proj2d = ptcl_native%get_rmat()
        do j = 1, BOX
            do i = 1, BOX
                proj2d(i,j,1) = proj2d(i,j,1) + sum(phantom(i,j,:)) + 0.2*real(i)/real(BOX) - 0.1*real(j)/real(BOX)
            end do
        end do
        call ptcl_native%set_rmat(proj2d, .false.)
        ! gridding route (prep_imgs4rec); obs_g retains the tapered crop
        call ptcl_work%copy(ptcl_native)
        call ptcl_work%norm_noise_fft_clip_shift(lmsk_native, obs_g, [0.,0.])
        call obs_g%ifft
        call obs_g%norm_noise_taper_edge_pad_fft(lmsk_crop, pad_g, renorm=.false.)
        allocate(obsg(CROP_BOX,CROP_BOX), obsp(CROP_BOX,CROP_BOX))
        rm_pad = obs_g%get_rmat()
        obsg   = rm_pad(:,:,1)
        ! PCG route
        call ptcl_work%copy(ptcl_native)
        call prep_rec_observation(ptcl_work, lmsk_native, obs_p, .true.)
        rm_pad = obs_p%get_rmat()
        obsp   = rm_pad(:,:,1)
        obs_scale = sum(obsg*obsp) / max(sum(obsp*obsp), TINY)
        obs_err   = sqrt(sum((obsg - obs_scale*obsp)**2) / max(sum(obsg**2), TINY))
        write(logfhandle,'(a,f10.6,a,es12.4)') '    gridding/PCG observation scale = ', obs_scale, &
            &'  scaled relative error = ', obs_err
        if( obs_err > OBS_PARITY_RELTOL )then
            write(logfhandle,'(a)') '    FAIL: gridding and PCG prepare different cropped observations'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: gridding and PCG prepare the same cropped observation'
        endif
        call ptcl_native%kill
        call ptcl_work%kill
        call obs_g%kill
        call obs_p%kill
        call pad_g%kill
        deallocate(proj2d, rm_pad, obsg, obsp, lmsk_native, lmsk_crop)
    else
        write(logfhandle,'(a)') '>>> STAGE 13 SKIPPED: an earlier stage failed'
    endif

    ! ============ STAGE 14: support semantics and warm-start band stability ============
    ! The shipped PCG map is window*u with u solved on the hard domain window > 0.
    ! Verify that the shipped map has no leakage outside that domain and that
    ! repeated output-space warm starts do not shrink the soft window band (an
    ! extra window multiplication on solver entry would).
    if( all_ok )then
        write(logfhandle,'(a)') '>>> STAGE 14: support semantics and warm-start band stability'
        call pcg_c%new(BOX, SMPD, LAMBDA)
        call pcg_c%set_deapod(.false.)
        call pcg_c%set_window_sphere(SUPPORT_MSKRAD)
        call pcg_c%prep_particles(projdirs, use_ctf=.true., sig2=sig2_2d)
        call pcg_c%begin_accum
        call pcg_c%accumulate_batch(y_planes, NPROJS, 1)
        call pcg_c%end_accum(.true.)
        call pcg_c%set_op_mode(PCG_OP_KERNEL)
        allocate(x_c(BOX,BOX,BOX), source=0.0)
        call pcg_c%solve_accum(x_c, maxits=SUPPORT_ITS, rtol=SUPPORT_RTOL, outcome=solver_outcome)
        write(logfhandle,'(a,a,a,i0,a,es12.4)') '    constrained solve: ', trim(solver_outcome%stop_reason), &
            &' after ', solver_outcome%iteration_count, ' iterations, residual = ', solver_outcome%final_rel_residual
        if( trim(solver_outcome%stop_reason) == PCG_STOP_INDEFINITE )then
            write(logfhandle,'(a)') '    FAIL: the constrained solve lost positive-definiteness'
            all_ok = .false.
        endif
        ! the window itself, by the set_window_sphere recipe
        allocate(window(BOX,BOX,BOX), source=1.0)
        call wimg%new([BOX,BOX,BOX], SMPD)
        call wimg%set_rmat(window, .false.)
        call wimg%mask3D_soft(SUPPORT_MSKRAD, backgr=0.)
        window = wimg%get_rmat()
        call wimg%kill
        support_leak = maxval(abs(x_c), mask=window <= TINY)
        write(logfhandle,'(a,es14.6)') '    maximum outside-support magnitude = ', support_leak
        if( support_leak <= TINY )then
            write(logfhandle,'(a)') '    PASS: the constrained output is zero outside its support'
        else
            write(logfhandle,'(a)') '    FAIL: the constrained output leaks outside its support'
            all_ok = .false.
        endif
        rms_first = -1.0
        rms_last  = -1.0
        do irep = 1, SUPPORT_NREP
            call pcg_c%solve_accum(x_c, maxits=2, rtol=0.0)
            call band_rms(x_c, rms_c)
            if( irep == 1 ) rms_first = rms_c
            rms_last = rms_c
            write(logfhandle,'(a,i3,a,es14.6)') '    warm start ', irep, ': band rms = ', rms_c
        end do
        if( rms_last < SUPPORT_STABILITY_FRAC * rms_first )then
            write(logfhandle,'(a)') '    FAIL: the window band shrinks under repeated warm starts'
            all_ok = .false.
        else
            write(logfhandle,'(a)') '    PASS: the window band is stable under repeated warm starts'
        endif
        call pcg_c%kill
        deallocate(x_c, window)
    else
        write(logfhandle,'(a)') '>>> STAGE 14 SKIPPED: an earlier stage failed'
    endif

    call pcgop%kill
    call pcg_reduce%kill
    call pcg_crop%kill
    call pcg_ml%kill
    call projdirs%kill
    call projdirs_exp%kill
    call projdirs_crop%kill
    call e%kill
    call e_exp%kill
    call c1sym%kill
    call c2sym%kill
    if( all_ok )then
        call simple_end('**** SIMPLE_TEST_PCG_RECON NORMAL STOP ****')
    else
        THROW_HARD('TEST_PCG_RECON FAILED')
    endif

  contains

    real function corr_of( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        real :: ma, mb, num, da, db
        ma  = sum(a)/real(size(a))
        mb  = sum(b)/real(size(b))
        num = sum((a-ma)*(b-mb))
        da  = sum((a-ma)**2)
        db  = sum((b-mb)**2)
        corr_of = num / sqrt(max(da*db, TINY))
    end function corr_of

    !> RMS of a volume over the window band [SUPPORT_BAND_LO, SUPPORT_BAND_HI]
    subroutine band_rms( vol, rms )
        real, intent(in)  :: vol(:,:,:)
        real, intent(out) :: rms
        integer :: n
        n = count(window >= SUPPORT_BAND_LO .and. window <= SUPPORT_BAND_HI)
        rms = 0.
        if( n > 0 ) rms = sqrt(sum(vol**2, mask=(window >= SUPPORT_BAND_LO .and. window <= SUPPORT_BAND_HI)) / real(n))
    end subroutine band_rms

end subroutine exec_test_pcg_recon

end submodule simple_commanders_test_highlevel_pcg
