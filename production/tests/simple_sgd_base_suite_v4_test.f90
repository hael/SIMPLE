module base_sgd_v4_test
    use base_sgd_test_helpers, only: prepare_test_workspace, make_1jyx_volume, &
                                     simulate_1jyx_particles, check_stack_count, create_import_project
contains
    subroutine run_abinitio_v4()
        use simple_core_module_api
        use simple_atoms, only: atoms
        use simple_molecule_data, only: molecule_data, betagal_1jyx
        use simple_cmdline, only: cmdline
        use simple_commanders_sim, only: commander_simulate_particles
        use simple_commanders_project_core, only: commander_new_project
        use simple_commanders_project_ptcl, only: commander_import_particles
        use simple_sp_project, only: sp_project
        use simple_parameters, only: parameters
        use simple_builder, only: builder
        use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
        use simple_type_defs, only: OBJFUN_EUCLID
        use simple_image, only: image
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
#include "simple_local_flags.inc"

        character(len=*), parameter :: WORKDIR = 'test_1jyx_abinitio_v4'
        character(len=*), parameter :: PROJNAME = 'onejyx_abinitio_v4'
        character(len=*), parameter :: PROJFILE = PROJNAME//'.simple'
        character(len=*), parameter :: VOLFILE = '1JYX_v4.mrc'
        character(len=*), parameter :: CLEANSTK = '1JYX_v4_clean.mrcs'
        character(len=*), parameter :: NOISYSTK = '1JYX_v4_noisy.mrcs'

        real :: smpd, mskdiam, snr, shift_limits(2, 2)
        integer :: nptcls, nthr, status, ldim(3), vol_dim(3)
        integer :: iref, irot, candidate_rot, best_ref, best_rot, truth_ref, truth_rot
        integer :: seed_size, iseed, candidate_count
        integer, allocatable :: seed(:)
        integer :: pdim_srch(3)
        real :: truth_angle, applied_shift(2), expected_shift(2), recovered(3)
        real :: angle_err, angle_err_alt, recovered_angle
        real(dp) :: loss, grad(2), best_loss, objective_initial, objective_final
        real(dp) :: candidate_initial, candidate_final
        real(dp) :: scalar0, scalar_xp, scalar_xm, scalar_fd_x, diag_delta, direct_score
        real(dp) :: vector0, vector_xp, vector_xm
        real(dp) :: diag_shifts(2,5), l_scalar_diag, l_vector_diag, l_direct_diag
        real(dp) :: dmin, dmax, dshift_min, dshift_max, direct_delta_max
        real(dp) :: vector_fd_x, scalar_grad_x, max_grad_error
        integer :: diag_rots(3), diag_ref, diag_rot_index, diag_shift_index
        integer :: scalar_best_rot, vector_best_rot
        real(dp) :: scalar_best_loss, vector_best_loss
        real :: candidate_result(3)
        integer :: accepted_steps, candidate_steps, best_steps
        logical :: pass_class, pass_angle, pass_loss, pass_shift, pass_roundtrip, pass_all
        real, allocatable, target :: sigma2_noise(:, :)
        complex(sp), allocatable :: ref_pft_diag(:, :), ptcl_pft_diag(:, :)
        real(sp), allocatable :: legacy_scores(:)
        real(sp), allocatable :: raw_vector_diag(:)
        type(cmdline) :: cpft
        type(commander_simulate_particles) :: xsim
        type(commander_new_project) :: xnew
        type(commander_import_particles) :: ximport
        type(molecule_data) :: mol
        type(atoms) :: molecule
        type(builder) :: b
        type(pftc_shsrch_grad) :: search
        type(parameters) :: params
        type(image) :: ref1, ref2, observed
        type(string) :: cwd, root, project_path, clean_path, noisy_path, vol_path

        smpd = 1.3; mskdiam = 120.; snr = 10.; nptcls = 2; nthr = 4; vol_dim = [144, 144, 144]
        ! Keep the truth-controlled fixture reproducible.  The same seed covers
        ! the simulated orientations and the subsequent Gaussian noise, so a
        ! failed candidate search is a repeatable numerical failure rather than
        ! an accidental draw from the global RNG state.
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        do iseed = 1, seed_size
            seed(iseed) = 20260730 + 97 * (iseed - 1)
        enddo
        call random_seed(put=seed)
        deallocate(seed)
        truth_ref = 2; truth_angle = 37.; applied_shift = [2., -1.5]
        ! rtsq applies the translation in the rotated image frame.  The production
        ! matcher returns the inverse (corrective) shift in that same frame, hence
        ! expected_shift = -R(truth_angle) applied_shift, not simply -applied_shift.
        expected_shift(1) = -(cos(deg2rad(truth_angle))*applied_shift(1) - sin(deg2rad(truth_angle))*applied_shift(2))
        expected_shift(2) = -(sin(deg2rad(truth_angle))*applied_shift(1) + cos(deg2rad(truth_angle))*applied_shift(2))
        call prepare_test_workspace(WORKDIR, cwd, root, status)

        write (logfhandle, '(a)') '>>> V4 STEP 1: generate volume and two distinct reference projections'
        mol = betagal_1jyx()
        call make_1jyx_volume(VOLFILE, smpd, mol, molecule, vol_path, vol_dim)
        call simulate_1jyx_particles(xsim, vol_path, CLEANSTK, smpd, mskdiam, nptcls, nthr, outfile='v4_simulated.simple')
        clean_path = simple_abspath(string(CLEANSTK))
        call check_stack_count(clean_path, nptcls, 'V4 clean stack', ldim)

        call ref1%new([ldim(1), ldim(2), 1], smpd); call ref1%read(clean_path, 1)
        call ref2%new([ldim(1), ldim(2), 1], smpd); call ref2%read(clean_path, 2)
        ! The synthetic particle is the second reference, rotated and translated by
        ! the known truth shift.  SIMPLE's alignment result is the corrective shift,
        ! so the expected recovered value is -applied_shift.
        call ref2%rtsq(truth_angle, applied_shift(1), applied_shift(2), observed); call observed%add_gauran(snr)
        call observed%write(string(NOISYSTK), 1, del_if_exists=.true.); noisy_path = simple_abspath(string(NOISYSTK))

        write (logfhandle, '(a)') '>>> V4 STEP 2: create one-particle project and production polar context'
        call create_import_project(xnew, ximport, PROJNAME, PROJFILE, noisy_path, smpd, project_path)
        call cpft%set('projfile', project_path); call cpft%set('smpd', smpd); call cpft%set('mskdiam', mskdiam); call cpft%set('ncls', 2); call cpft%set('nptcls', 1); call cpft%set('ctf', 'no'); call cpft%set('lp', 8.)
        call b%init_params_and_build_strategy2D_tbox(cpft, params, wthreads=.false.); params%cc_objfun = OBJFUN_EUCLID
        call b%pftc%new(params, 2, [1, 1], params%kfromto)
        call ref1%fft(); call ref2%fft(); call observed%fft()
        call ref1%memoize4polarize(b%pftc%get_pdim_srch()); call ref2%memoize4polarize(b%pftc%get_pdim_srch()); call observed%memoize4polarize(b%pftc%get_pdim_srch())
        call b%pftc%polarize_ref_pft(ref1, 1, .true., b%pftc%get_pdim_srch(), .false.); call b%pftc%polarize_ref_pft(ref2, 2, .true., b%pftc%get_pdim_srch(), .false.)
        call b%pftc%polarize_ptcl_pft(observed, 1, b%pftc%get_pdim_srch(), .false.); call b%pftc%set_eo(1, .true.)
        ! The scalar gradient consumes pfts_refs_* directly, whereas the
        ! vector raw-loss path consumes the memoized FT(ref) arrays.  Populate
        ! both representations before the round-trip diagnostic and candidate
        ! search; otherwise the zero-shift vector call sees uninitialized data.
        call b%pftc%memoize_refs
        allocate (sigma2_noise(params%kfromto(1):params%kfromto(2), 1), source=0.05); call b%pftc%assign_sigma2_noise(sigma2_noise)
        ! The vector path also consumes the memoized particle/CTF products
        ! (ft_ptcl_ctf and ft_ctf2); the scalar gradient computes these terms
        ! directly, so populate this second production memo before comparing
        ! the two APIs.
        call b%pftc%memoize_ptcls
        ! polarize_ptcl_pft memoizes the weighted norm at polarization time.  Sigma
        ! calibration is attached immediately afterward here, so refresh that cache
        ! before evaluating the raw Euclidean loss; otherwise its denominator is zero.
        call b%pftc%memoize_sqsum_ptcl(1)
        pdim_srch = b%pftc%get_pdim_srch()
        allocate (ref_pft_diag(pdim_srch(1), pdim_srch(2):pdim_srch(3)))
        allocate (ptcl_pft_diag(pdim_srch(1), pdim_srch(2):pdim_srch(3)))
        call b%pftc%get_ref_pft(1, .true., ref_pft_diag)
        write (logfhandle, '(a,es16.8)') '>>> V4 POLAR ENERGY REF1: ', sum(real(ref_pft_diag*conjg(ref_pft_diag), dp))
        call b%pftc%get_ref_pft(2, .true., ref_pft_diag)
        write (logfhandle, '(a,es16.8)') '>>> V4 POLAR ENERGY REF2: ', sum(real(ref_pft_diag*conjg(ref_pft_diag), dp))
        call b%pftc%get_ptcl_pft(1, ptcl_pft_diag)
        write (logfhandle, '(a,es16.8)') '>>> V4 POLAR ENERGY PTCL: ', sum(real(ptcl_pft_diag*conjg(ptcl_pft_diag), dp))
        deallocate (ref_pft_diag, ptcl_pft_diag)
        ! Match the production calibration boundary: the raw Euclidean objective
        ! requires a finite per-shell variance, not merely an allocated array.
        call b%pftc%gen_sigma_contrib(2, 1, [0.0, 0.0], 1, sigma2_noise(:, 1))
        sigma2_noise(:, 1) = max(sigma2_noise(:, 1), 1.0e-6)
        write (logfhandle, '(a,2es16.8)') '>>> V4 SIGMA2 MIN/MAX: ', minval(sigma2_noise(:, 1)), maxval(sigma2_noise(:, 1))
        ! Probe each reference before searching the full rotation grid.  These values
        ! are the raw finite loss L and its shift gradient at s=(0,0); NaN here means
        ! the production Fourier/polar context is invalid before SGD is entered.
        call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0._dp, 0._dp], 1, loss, grad)
        write (logfhandle, '(a,2es16.8,1x,l1)') '>>> V4 RAW PROBE REF1 (LOSS,GX): ', loss, grad(1), ieee_is_finite(loss)
        ! Diagnostic phase only: evaluate the scalar direct-gradient objective
        ! at a symmetric finite-difference stencil and compare it with the
        ! established FFT/vector raw-loss API at the same states.  This does
        ! not alter optimization; it distinguishes a gradient/sign problem
        ! from a residual normalization or rotation-indexing problem.
        ! Use a displacement larger than sqrt(SHERRSQ) (about 0.0032 px in
        ! SIMPLE).  Below that threshold the vector implementation is allowed
        ! to reuse its zero-shift memoized FFT state, which would make the
        ! finite-difference probe appear falsely constant.
        diag_delta = 1.0e-2_dp
        scalar0 = loss
        call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [diag_delta, 0._dp], 1, scalar_xp, grad)
        call b%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [-diag_delta, 0._dp], 1, scalar_xm, grad)
        scalar_fd_x = (scalar_xp - scalar_xm) / (2._dp * diag_delta)
        allocate(raw_vector_diag(b%pftc%get_nrots()))
        call b%pftc%gen_raw_euclid_vals(1, 1, [0._sp, 0._sp], raw_vector_diag)
        vector0 = real(raw_vector_diag(1), dp)
        call b%pftc%gen_raw_euclid_vals(1, 1, [real(diag_delta,sp), 0._sp], raw_vector_diag)
        vector_xp = real(raw_vector_diag(1), dp)
        call b%pftc%gen_raw_euclid_vals(1, 1, [-real(diag_delta,sp), 0._sp], raw_vector_diag)
        write(logfhandle,'(a,4es18.10)') '>>> V4 DIAG SCALAR 0/+X/-X/FDX: ', &
            scalar0, scalar_xp, scalar_xm, scalar_fd_x
        vector_xm = real(raw_vector_diag(1), dp)
        write(logfhandle,'(a,3es18.10)') '>>> V4 DIAG VECTOR 0/+X/-X: ', &
            vector0, vector_xp, vector_xm
        direct_score = b%pftc%gen_corr_for_rot_8(1, 1, [0._dp, 0._dp], 1)
        if (ieee_is_finite(direct_score) .and. direct_score > 0.0_dp) then
        write(logfhandle,'(a,2es18.10)') '>>> V4 DIAG DIRECT SCORE/RAW ROT1: ', &
                direct_score, -log(direct_score)
        else
            write(logfhandle,'(a,es18.10)') '>>> V4 DIAG DIRECT SCORE ROT1: ', direct_score
        endif
        deallocate(raw_vector_diag)

        ! Characterize the objective difference over a small deterministic
        ! state grid.  This is diagnostic only: it does not select candidates
        ! or modify the optimizer.  The shift step is above SHERRSQ so the
        ! vector path recomputes its shifted FFT state rather than reusing the
        ! zero-shift memo.
        diag_rots = [1, max(1, b%pftc%get_nrots()/2), b%pftc%get_nrots()]
        diag_shifts(:,1) = [0._dp, 0._dp]
        diag_shifts(:,2) = [0.01_dp, 0._dp]
        diag_shifts(:,3) = [-0.01_dp, 0._dp]
        diag_shifts(:,4) = [0._dp, 0.01_dp]
        diag_shifts(:,5) = [0._dp, -0.01_dp]
        dmin = huge(1._dp); dmax = -huge(1._dp)
        dshift_min = huge(1._dp); dshift_max = -huge(1._dp)
        direct_delta_max = 0._dp; max_grad_error = 0._dp
        allocate(raw_vector_diag(b%pftc%get_nrots()))
        do diag_ref = 1, 2
            do diag_rot_index = 1, size(diag_rots)
                dshift_min = huge(1._dp); dshift_max = -huge(1._dp)
                do diag_shift_index = 1, size(diag_shifts,2)
                    call b%pftc%gen_raw_euclid_grad_for_rot_8(diag_ref, 1, &
                        diag_shifts(:,diag_shift_index), diag_rots(diag_rot_index), l_scalar_diag, grad)
                    call b%pftc%gen_raw_euclid_vals(diag_ref, 1, &
                        real(diag_shifts(:,diag_shift_index),sp), raw_vector_diag)
                    l_vector_diag = real(raw_vector_diag(diag_rots(diag_rot_index)),dp)
                    direct_score = b%pftc%gen_corr_for_rot_8(diag_ref, 1, &
                        diag_shifts(:,diag_shift_index), diag_rots(diag_rot_index))
                    if (direct_score > 0._dp .and. ieee_is_finite(direct_score)) then
                        l_direct_diag = -log(direct_score)
                    else
                        l_direct_diag = huge(1._dp)
                    endif
                    dmin = min(dmin, l_vector_diag-l_scalar_diag)
                    dmax = max(dmax, l_vector_diag-l_scalar_diag)
                    dshift_min = min(dshift_min, l_vector_diag-l_scalar_diag)
                    dshift_max = max(dshift_max, l_vector_diag-l_scalar_diag)
                    direct_delta_max = max(direct_delta_max, abs(l_direct_diag-l_scalar_diag))
                    if (diag_ref == 1 .and. diag_rot_index == 1 .and. diag_shift_index == 1) then
                        scalar_grad_x = grad(1)
                        call b%pftc%gen_raw_euclid_vals(diag_ref, 1, [real(diag_delta,sp),0._sp], raw_vector_diag)
                        vector_xp = real(raw_vector_diag(diag_rots(diag_rot_index)),dp)
                        call b%pftc%gen_raw_euclid_vals(diag_ref, 1, [-real(diag_delta,sp),0._sp], raw_vector_diag)
                        vector_xm = real(raw_vector_diag(diag_rots(diag_rot_index)),dp)
                        vector_fd_x = (vector_xp-vector_xm)/(2._dp*diag_delta)
                        max_grad_error = abs(vector_fd_x-scalar_grad_x)
                    endif
                enddo
                write(logfhandle,'(a,2i5,2es18.10)') '>>> V4 DIAG CANDIDATE (REF,ROT,DSHIFT_MIN,DSHIFT_MAX): ', &
                    diag_ref, diag_rots(diag_rot_index), dshift_min, dshift_max
            enddo
        enddo
        deallocate(raw_vector_diag)
        write(logfhandle,'(a,3es18.10)') '>>> V4 DIAG VECTOR-SCALAR D MIN/MAX/RANGE: ', &
            dmin, dmax, dmax-dmin
        write(logfhandle,'(a,2es18.10)') '>>> V4 DIAG DIRECT-SCALAR MAX/DX: ', &
            direct_delta_max, max_grad_error
        ! Compare the discrete zero-shift winner under both raw APIs.
        allocate(raw_vector_diag(b%pftc%get_nrots()))
        do diag_ref = 1, 2
            scalar_best_loss = huge(1._dp); vector_best_loss = huge(1._dp)
            scalar_best_rot = 0; vector_best_rot = 0
            do diag_rot_index = 1, b%pftc%get_nrots()
                call b%pftc%gen_raw_euclid_grad_for_rot_8(diag_ref, 1, [0._dp,0._dp], &
                    diag_rot_index, l_scalar_diag, grad)
                call b%pftc%gen_raw_euclid_vals(diag_ref, 1, [0._sp,0._sp], raw_vector_diag)
                l_vector_diag = real(raw_vector_diag(diag_rot_index),dp)
                if (l_scalar_diag < scalar_best_loss) then
                    scalar_best_loss = l_scalar_diag; scalar_best_rot = diag_rot_index
                endif
                if (l_vector_diag < vector_best_loss) then
                    vector_best_loss = l_vector_diag; vector_best_rot = diag_rot_index
                endif
            enddo
            write(logfhandle,'(a,i4,2(i6,1x),2es18.10)') '>>> V4 DIAG ZERO-SHIFT ARGMIN REF/SCALAR/VECTOR: ', &
                diag_ref, scalar_best_rot, vector_best_rot, scalar_best_loss, vector_best_loss
        enddo
        deallocate(raw_vector_diag)
        call b%pftc%gen_raw_euclid_grad_for_rot_8(2, 1, [0._dp, 0._dp], 1, loss, grad)
        write (logfhandle, '(a,2es16.8,1x,l1)') '>>> V4 RAW PROBE REF2 (LOSS,GX): ', loss, grad(1), ieee_is_finite(loss)
        allocate (legacy_scores(b%pftc%get_nrots()))
        call b%pftc%gen_objfun_vals(1, 1, [0.0, 0.0], legacy_scores)
        if (legacy_scores(1) > 0.0 .and. ieee_is_finite(legacy_scores(1))) then
            write (logfhandle, '(a,2es16.8)') '>>> V4 LEGACY SCORE/RAW REF1 ROT1: ', legacy_scores(1), -log(real(legacy_scores(1), dp))
        else
            write (logfhandle, '(a,es16.8)') '>>> V4 LEGACY SCORE REF1 ROT1: ', legacy_scores(1)
        endif
        deallocate (legacy_scores)

        write (logfhandle, '(a)') '>>> V4 STEP 3: joint class/rotation search with bounded shift SGD'
        best_loss = huge(1._dp); best_ref = 0; best_rot = 0; best_steps = 0; candidate_count = 0
        write (logfhandle, '(a,i0)') '>>> V4 ROTATION COUNT: ', b%pftc%get_nrots()
        shift_limits(:, 1) = -5.; shift_limits(:, 2) = 5.
        call search%new_direct(b, shift_limits)
        call search%set_diagnostic_mode(.true.)
        do iref = 1, 2; do irot = 1, b%pftc%get_nrots()
                call b%pftc%gen_raw_euclid_grad_for_rot_8(iref, 1, [0._dp, 0._dp], irot, loss, grad)
                if (ieee_is_finite(loss)) then
                    candidate_count = candidate_count + 1
                    call search%set_indices(iref, 1)
                    candidate_rot = irot
                    candidate_result = search%minimize_direct(candidate_rot, [0.0, 0.0], .5, 8, sh_rot=.false., &
                        accepted_steps=candidate_steps, objective_initial=candidate_initial, &
                        objective_final=candidate_final, raw_euclid=.true.)
                    if (ieee_is_finite(candidate_final) .and. candidate_final < best_loss) then
                        best_loss = candidate_final
                        best_ref = iref
                        best_rot = candidate_rot
                        recovered = candidate_result
                        objective_initial = candidate_initial
                        objective_final = candidate_final
                        best_steps = candidate_steps
                    end if
                end if
            end do; end do
        if (best_rot < 1) THROW_HARD('V4 discrete class/rotation search produced no finite candidate')
        accepted_steps = best_steps
        write (logfhandle, '(a,i0,a,es16.8)') '>>> V4 JOINT CANDIDATES: ', candidate_count, ' BEST FINAL LOSS: ', best_loss
        ! Alignment reports the corrective rotation, so compare against -truth_angle.
        truth_rot = b%pftc%get_roind(real(-truth_angle, sp))
        call search%set_indices(best_ref, 1)
        irot = best_rot
        write (logfhandle, '(a,2i8)') '>>> V4 CLASS TRUE/RECOVERED: ', truth_ref, best_ref
        write (logfhandle, '(a,2i8)') '>>> V4 ROTATION TRUE/RECOVERED: ', b%pftc%get_roind(real(truth_angle, sp)), irot
        write (logfhandle, '(a,2i8)') '>>> V4 ROTATION CORRECTIVE/RECOVERED: ', truth_rot, irot
        write (logfhandle, '(a,2f10.4)') '>>> V4 ANGLE TRUE/RECOVERED: ', truth_angle, b%pftc%get_rot(irot)
        write (logfhandle, '(a,"(",f10.4,1x,f10.4,") (",f10.4,1x,f10.4,")")') &
            '>>> V4 SHIFT EXPECTED/RECOVERED: ', expected_shift(1), expected_shift(2), recovered(2), recovered(3)
        recovered_angle = b%pftc%get_rot(irot)
        angle_err = abs(truth_angle - recovered_angle); if (angle_err > 180.) angle_err = 360.-angle_err
        ! SIMPLE's in-plane rotation convention may report the corrective angle,
        ! i.e. the negative of the synthetic angle.  Circularize both alternatives
        ! before selecting the closer one.
        angle_err_alt = abs(-truth_angle - recovered_angle); if (angle_err_alt > 180.) angle_err_alt = 360.-angle_err_alt
        angle_err = min(angle_err, angle_err_alt)
        pass_class = best_ref == truth_ref
        pass_angle = angle_err <= max(2.*b%pftc%get_dang(), 5.)
        pass_loss = accepted_steps >= 1 .and. objective_final < objective_initial
        pass_shift = maxval(abs(real(recovered(2:3), dp) - real(expected_shift, dp))) <= 0.5
        write (logfhandle, '(a,l1)') '>>> V4 CHECK CLASS: ', pass_class
        write (logfhandle, '(a,l1)') '>>> V4 CHECK ANGLE: ', pass_angle
        write (logfhandle, '(a,l1)') '>>> V4 CHECK LOSS:  ', pass_loss
        write (logfhandle, '(a,l1)') '>>> V4 CHECK SHIFT: ', pass_shift
        pass_roundtrip = .not. search%diagnostic_failed()
        write (logfhandle, '(a,l1)') '>>> V4 CHECK ROUNDTRIP: ', pass_roundtrip
        pass_all = pass_class .and. pass_angle .and. pass_loss .and. pass_shift .and. pass_roundtrip
        call search%kill; call b%kill_strategy2D_tbox; call b%kill_general_tbox; deallocate (sigma2_noise)
        call simple_chdir(cwd, status); if (status /= 0) THROW_HARD('could not restore original working directory')
        if (.not. pass_all) THROW_HARD('V4 truth-controlled assignment regression failed; see CHECK lines')
        write (logfhandle, '(a)') '>>> V4 RESULTS: '//root%to_char(); call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO_V4 NORMAL STOP ****')
    end subroutine run_abinitio_v4
end module base_sgd_v4_test
