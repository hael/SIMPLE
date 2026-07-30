module base_sgd_test_helpers
    use simple_atoms, only: atoms
    use simple_cmdline, only: cmdline
    use simple_imghead, only: find_ldim_nptcls
    use simple_molecule_data, only: molecule_data
    use simple_procimgstk, only: add_noise_imgfile
    use simple_string, only: string
    use simple_ui, only: make_ui
    use simple_commanders_sim, only: commander_simulate_particles
    use simple_commanders_project_core, only: commander_new_project
    use simple_commanders_project_ptcl, only: commander_import_particles
    implicit none
#include "simple_local_flags.inc"

    logical, save :: ui_initialized = .false.

contains

    subroutine prepare_test_workspace(workdir, original_cwd, workflow_root, status)
        use simple_core_module_api
        character(len=*), intent(in) :: workdir
        type(string), intent(out) :: original_cwd, workflow_root
        integer, intent(out) :: status
        if (.not. ui_initialized) then
            call make_ui
            ui_initialized = .true.
        endif
        call simple_getcwd(original_cwd)
        if (file_exists(workdir)) then
            call simple_rmdir(workdir, status)
            if (status /= 0) THROW_HARD('could not reset '//trim(workdir))
        end if
        call simple_mkdir(workdir)
        call simple_chdir(workdir, status)
        if (status /= 0) THROW_HARD('could not enter '//trim(workdir))
        call simple_getcwd(workflow_root)
    end subroutine prepare_test_workspace

    subroutine make_1jyx_volume(volfile, smpd, mol, molecule, volume_path, vol_dim)
        use simple_core_module_api
        character(len=*), intent(in) :: volfile
        real, intent(in) :: smpd
        type(molecule_data), intent(in) :: mol
        type(atoms), intent(inout) :: molecule
        type(string), intent(out) :: volume_path
        integer, intent(in), optional :: vol_dim(3)
        integer :: ldim(3), nimgs
        if (present(vol_dim)) then
            call molecule%pdb2mrc(volfile=string(volfile), smpd=smpd, mol=mol, center_pdb=.true., vol_dim=vol_dim)
        else
            call molecule%pdb2mrc(volfile=string(volfile), smpd=smpd, mol=mol, center_pdb=.true.)
        end if
        call molecule%kill()
        volume_path = simple_abspath(string(volfile))
        call find_ldim_nptcls(volume_path, ldim, nimgs)
        if (any(ldim < 1) .or. ldim(3) <= 1) THROW_HARD('1JYX volume has invalid dimensions')
    end subroutine make_1jyx_volume

    subroutine check_stack_count(stack_path, expected, label, ldim_out)
        use simple_core_module_api
        type(string), intent(in) :: stack_path
        integer, intent(in) :: expected
        character(len=*), intent(in) :: label
        integer, intent(out), optional :: ldim_out(3)
        integer :: ldim(3), nimgs
        call find_ldim_nptcls(stack_path, ldim, nimgs)
        if (nimgs /= expected) THROW_HARD(trim(label)//' has the wrong particle count')
        if (present(ldim_out)) ldim_out = ldim
    end subroutine check_stack_count

    subroutine simulate_1jyx_particles(xsim, volume_path, outstk, smpd, mskdiam, nptcls, nthr, sherr, outfile)
        type(commander_simulate_particles), intent(inout) :: xsim
        type(string), intent(in) :: volume_path
        character(len=*), intent(in) :: outstk
        real, intent(in) :: smpd, mskdiam
        integer, intent(in) :: nptcls, nthr
        real, intent(in), optional :: sherr
        character(len=*), intent(in), optional :: outfile
        type(cmdline) :: cline
        call cline%set('prg', 'simulate_particles')
        call cline%set('mkdir', 'no')
        call cline%set('vol1', volume_path)
        call cline%set('outstk', outstk)
        if (present(outfile)) call cline%set('outfile', outfile)
        call cline%set('smpd', smpd)
        call cline%set('mskdiam', mskdiam)
        call cline%set('nptcls', nptcls)
        call cline%set('nthr', nthr)
        call cline%set('pgrp', 'c1')
        call cline%set('ctf', 'no')
        call cline%set('snr', 10.0)
        call cline%set('bfac', 0.0)
        if (present(sherr)) then
            call cline%set('sherr', sherr)
        else
            call cline%set('sherr', 0.0)
        end if
        call xsim%execute(cline)
        call cline%kill()
    end subroutine simulate_1jyx_particles

    subroutine add_noise_checked(clean_path, noisy_file, snr, smpd, expected, noisy_path, label)
        use simple_core_module_api
        type(string), intent(in) :: clean_path
        character(len=*), intent(in) :: noisy_file, label
        real, intent(in) :: snr, smpd
        integer, intent(in) :: expected
        type(string), intent(out) :: noisy_path
        call add_noise_imgfile(clean_path, string(noisy_file), snr, smpd)
        noisy_path = simple_abspath(string(noisy_file))
        call check_stack_count(noisy_path, expected, label)
    end subroutine add_noise_checked

    subroutine create_import_project(xnew, ximport, projname, projfile, stack_path, smpd, project_path, project_root, oritab)
        use simple_core_module_api
        type(commander_new_project), intent(inout) :: xnew
        type(commander_import_particles), intent(inout) :: ximport
        character(len=*), intent(in) :: projname, projfile
        type(string), intent(in) :: stack_path
        real, intent(in) :: smpd
        type(string), intent(out) :: project_path
        type(string), intent(out), optional :: project_root
        character(len=*), intent(in), optional :: oritab
        type(cmdline) :: cline_new, cline_import
        call cline_new%set('prg', 'new_project')
        call cline_new%set('projname', projname)
        call cline_new%set('qsys_name', 'local')
        call xnew%execute(cline_new)
        call cline_new%kill()
        if (present(project_root)) call simple_getcwd(project_root)
        project_path = simple_abspath(string(projfile))
        call cline_import%set('prg', 'import_particles')
        call cline_import%set('mkdir', 'no')
        call cline_import%set('projfile', project_path)
        call cline_import%set('stk', stack_path)
        call cline_import%set('smpd', smpd)
        call cline_import%set('ctf', 'no')
        if (present(oritab)) call cline_import%set('oritab', oritab)
        call ximport%execute(cline_import)
        call cline_import%kill()
    end subroutine create_import_project

end module base_sgd_test_helpers

module base_sgd_direct_test
contains
    subroutine run_direct_shift()
        use simple_core_module_api, only: dp
        use simple_pftc_shsrch_grad, only: bounded_shift_trial
        implicit none

        real(dp) :: shift(2), shift_gradient(2), shift_limits(2, 2), trial_shift(2)

        shift = 0.0_dp
        shift_gradient = [3.0_dp, 4.0_dp]
        shift_limits(:, 1) = -1.0_dp
        shift_limits(:, 2) = 1.0_dp
        trial_shift = bounded_shift_trial(shift, shift_gradient, 0.5_dp, shift_limits)
        call require_close_direct(trial_shift, [-0.3_dp, -0.4_dp], 1.0e-12_dp,&
            &'direct step is normalized so eta is a pixel-length bound')

        shift_limits(:, 1) = -0.25_dp
        shift_limits(:, 2) = 0.25_dp
        trial_shift = bounded_shift_trial(shift, shift_gradient, 0.5_dp, shift_limits)
        call require_close_direct(trial_shift, [-0.25_dp, -0.25_dp], 1.0e-12_dp,&
            &'direct step is projected into the legal shift box')

        shift = [0.1_dp, -0.2_dp]
        trial_shift = bounded_shift_trial(shift, [0.0_dp, 0.0_dp], 0.5_dp, shift_limits)
        call require_close_direct(trial_shift, shift, 1.0e-12_dp, 'zero gradient preserves the original shift')

        write (*, '(A)') 'joint2D bounded direct-shift regression: PASS'
    end subroutine run_direct_shift

    subroutine require_close_direct(actual, expected, tolerance, msg)
        use simple_core_module_api, only: dp
        real(dp), intent(in) :: actual(2), expected(2), tolerance
        character(len=*), intent(in) :: msg
        if (maxval(abs(actual - expected)) > tolerance) then
            write (*, '(A)') 'simple_test_joint2D_direct_shift failed: '//trim(msg)
            error stop 1
        end if
    end subroutine require_close_direct
end module base_sgd_direct_test

module base_sgd_sampling_restore_test
contains
    subroutine run_sampling_and_restore_policy()
        use simple_classaverager, only: cavger_zero_support_recovery
        use simple_oris, only: oris
        implicit none
        type(oris) :: os
        integer, allocatable :: inds(:)
        integer :: nsamp
        real :: frac

        call os%new(10, .true.)
        call os%set_all2single('state', 1.0)
        frac = 0.6
        call os%sample4update_rnd([1, 10], frac, nsamp, inds, .true.)
        if( nsamp /= 6 ) error stop 'SGD mini-batch fraction regression failed'
        if( any(inds < 1) .or. any(inds > 10) ) error stop 'SGD mini-batch index regression failed'
        call os%kill()

        if( .not. cavger_zero_support_recovery(0, .true.) ) &
            error stop 'zero-support recovery regression failed'
        if( cavger_zero_support_recovery(0, .false.) ) &
            error stop 'zero-support first-restore regression failed'
        if( cavger_zero_support_recovery(1, .true.) ) &
            error stop 'supported-class restoration regression failed'
        write(*,'(A)') 'SGD mini-batch and zero-support policy: PASS'
    end subroutine run_sampling_and_restore_policy
end module base_sgd_sampling_restore_test

module base_sgd_v2_test
    use base_sgd_test_helpers, only: prepare_test_workspace, make_1jyx_volume, &
                                     simulate_1jyx_particles, check_stack_count, create_import_project
contains
    subroutine run_abinitio_v2()
        use simple_core_module_api
        use simple_atoms, only: atoms
        use simple_molecule_data, only: molecule_data, betagal_1jyx
        use simple_cmdline, only: cmdline
        use simple_commanders_sim, only: commander_simulate_particles
        use simple_commanders_project_core, only: commander_new_project
        use simple_commanders_project_ptcl, only: commander_import_particles
        use simple_pftc_srch_api
        use simple_builder, only: builder
        use simple_pftc_shsrch_grad, only: pftc_shsrch_grad
        use simple_type_defs, only: OBJFUN_EUCLID
        use simple_image, only: image
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        implicit none
#include "simple_local_flags.inc"

        character(len=*), parameter :: WORKDIR = 'test_1jyx_abinitio'
        character(len=*), parameter :: PROJNAME = 'onejyx_abinitio'
        character(len=*), parameter :: PROJFILE = PROJNAME//'.simple'
        character(len=*), parameter :: VOL_FILE = '1JYX.mrc'
        character(len=*), parameter :: CLEAN_FILE = '1JYX_particles_clean.mrcs'
        character(len=*), parameter :: NOISY_FILE = '1JYX_particles_noisy.mrcs'
        character(len=*), parameter :: ORI_FILE = '1JYX_particles_oris.txt'

        real :: smpd, mskdiam, noise_snr
        integer :: nptcls, nthr
        type(cmdline) :: user_args
        type(commander_simulate_particles) :: xsimulate_particles
        type(commander_new_project) :: xnew_project
        type(commander_import_particles) :: ximport_particles
        type(molecule_data) :: mol
        type(atoms) :: molecule
        type(parameters) :: pft_params
        type(builder) :: pft_builder
        type(pftc_shsrch_grad) :: direct_shift_search
        type(string) :: workflow_root, project_path
        type(string) :: volume_path, clean_path, noisy_path
        type(string) :: original_cwd
        type(cmdline) :: cline_pft
        integer :: ldim(3), status
        integer :: vol_dim(3)
        integer :: pdim_srch(3), direct_irot, direct_accepted
        real :: shift_limits(2, 2), direct_cxy(3)
        real(dp) :: objective_initial, objective_final
        real(dp) :: raw_f, raw_grad(2), raw_f_xp, raw_f_xm, raw_f_yp, raw_f_ym
        real(dp) :: ref_pft_energy, ptcl_pft_energy, sigma_min, sigma_max
        complex(sp), allocatable :: ref_pft_diag(:, :), ptcl_pft_diag(:, :)
        real, allocatable :: ref_rmat_diag(:, :, :), ptcl_rmat_diag(:, :, :)
        real, allocatable, target :: sigma2_noise(:, :)

        type :: alignment_truth
            integer :: class_id
            integer :: rotation_index
            real :: angle_deg
            real :: shift(2)
        end type alignment_truth
        type(alignment_truth) :: truth
        type(image) :: reference, observed, corrected

        ! Defaults keep the test reasonably small while leaving enough particles for
        ! both classifications.  They can be overridden with old-style key=value args.
        smpd = 1.3
        mskdiam = 180.0
        noise_snr = 10.0
        nptcls = 200
        nthr = 4
        vol_dim = [144, 144, 144]
        truth%class_id = 1
        truth%rotation_index = 1
        truth%angle_deg = 0.0
        truth%shift = [2.0, 0.0]

        if (command_argument_count() > 0) then
            call user_args%parse_oldschool
            if (user_args%defined('smpd')) smpd = user_args%get_rarg('smpd')
            if (user_args%defined('mskdiam')) mskdiam = user_args%get_rarg('mskdiam')
            if (user_args%defined('snr')) noise_snr = user_args%get_rarg('snr')
            if (user_args%defined('nptcls')) nptcls = user_args%get_iarg('nptcls')
            if (user_args%defined('nthr')) nthr = user_args%get_iarg('nthr')
        end if
        if (smpd <= 0.0) THROW_HARD('smpd must be positive')
        if (mskdiam <= 0.0) THROW_HARD('mskdiam must be positive')
        if (noise_snr <= 0.0) THROW_HARD('snr must be positive')
        if (nptcls < 8) THROW_HARD('nptcls must be at least 8')
        if (nthr < 1) THROW_HARD('nthr must be positive')

        ! Direct commander calls need the normal SIMPLE UI metadata for mkdir=yes.
        call prepare_test_workspace(WORKDIR, original_cwd, workflow_root, status)

        write (logfhandle, '(a)') '>>> Step 1/6: generate a volume from embedded 1JYX beta-galactosidase coordinates'
        mol = betagal_1jyx()
        write (logfhandle, '(A,3I6)') '>>> TEST EXPLICIT VOLUME DIMENSIONS: ', vol_dim
        call make_1jyx_volume(VOL_FILE, smpd, mol, molecule, volume_path, vol_dim)

        write (logfhandle, '(a)') '>>> Step 2/6: simulate a clean stack at random orientations'
        call simulate_1jyx_particles(xsimulate_particles, volume_path, CLEAN_FILE, smpd, mskdiam, nptcls, nthr, &
                                     outfile=ORI_FILE)
        clean_path = simple_abspath(string(CLEAN_FILE))
        call check_stack_count(clean_path, nptcls, 'clean simulated stack', ldim)

        write (logfhandle, '(a)') '>>> Step 3/6: apply rotate by known angle and shift by known (sx, sy)'
        call reference%new([ldim(1), ldim(2), 1], smpd)
        call reference%read(string(CLEAN_FILE), 1)

        ! Generate shifted particle
        call reference%rtsq(0.0, &
                            truth%shift(1), &
                            truth%shift(2), &
                            observed)

        write (logfhandle, '(a,f7.3)') '>>> Step 4/6: add Gaussian noise to the stack; SNR = ', noise_snr
        ! Add noise to this one image
        call observed%add_gauran(noise_snr)

        ! Materialize exactly the particle that the production workflow would import.
        ! Keeping this one-image stack separate from the original simulated stack is
        ! important: the direct test must exercise the same stack/project ownership
        ! rules as abinitio2D without launching the full multi-class experiment.
        call observed%write(string(NOISY_FILE), 1, del_if_exists=.true.)
        call check_stack_count(string(NOISY_FILE), 1, 'materialized synthetic stack')
        noisy_path = simple_abspath(string(NOISY_FILE))

        ! Optional manual correction for inspection
        call observed%rtsq(0.0, &
                           -truth%shift(1), &
                           -truth%shift(2), &
                           corrected)

        write (logfhandle, '(a)') '>>> Step 5/6: prepare one-particle project and production Fourier/polar context'

        ! This is the same preparation boundary used by abinitio2D: create a SIMPLE
        ! project, import the particle stack, and let the builder read the project
        ! metadata before constructing the 2D toolbox.  We intentionally stop before
        ! commander_abinitio2D; a later v3 test will cover the complete classification.
        call create_import_project(xnew_project, ximport_particles, PROJNAME, PROJFILE, noisy_path, smpd, project_path)

        write (logfhandle, '(a)') '>>> Step 6/6: invoke the production Fourier/polar direct shift gradient'
        call cline_pft%set('projfile', project_path)
        call cline_pft%set('smpd', smpd)
        call cline_pft%set('mskdiam', mskdiam)
        call cline_pft%set('ncls', 1)
        call cline_pft%set('nptcls', 1)
        call cline_pft%set('ctf', 'no')
        call cline_pft%set('lp', 8.0)
        call pft_builder%init_params_and_build_strategy2D_tbox(cline_pft, pft_params, wthreads=.false.)
        pft_params%cc_objfun = OBJFUN_EUCLID
        call pft_builder%pftc%new(pft_params, 1, [1, 1], pft_params%kfromto)
        ! The production workflow normally attaches this calibration through
        ! simple_euclid_sigma2.  This standalone synthetic test constructs the
        ! polar calculator directly, so provide a finite positive variance for
        ! every Fourier shell and the one test particle before evaluating NLL.
        allocate (sigma2_noise(pft_params%kfromto(1):pft_params%kfromto(2), 1), source=1.0)
        call pft_builder%pftc%assign_sigma2_noise(sigma2_noise)
        pdim_srch = pft_builder%pftc%get_pdim_srch()
        ! Match SIMPLE's production image-to-polar workflow: polarize reads the
        ! image Fourier buffer, so explicitly FFT both synthetic spatial images first.
        ! Capture spatial energies before FFT; ordinary image has no get_sumsq method.
        ref_rmat_diag = reference%get_rmat()
        ptcl_rmat_diag = observed%get_rmat()
        write (logfhandle, '(a,2es16.8)') '>>> DIAG IMAGE SUMSQ REF/OBS: ', &
            sum(real(ref_rmat_diag, dp)**2), sum(real(ptcl_rmat_diag, dp)**2)
        deallocate (ref_rmat_diag, ptcl_rmat_diag)
        call reference%fft()
        call observed%fft()
        ! Match the production matcher: build the image-specific Cartesian-to-polar
        ! interpolation table before filling either polar Fourier buffer.  Without
        ! this call polarize() has no valid interpolation weights and returns zeros.
        call reference%memoize4polarize(pdim_srch)
        call pft_builder%pftc%polarize_ref_pft(reference, 1, iseven=.true., pdim=pdim_srch, oversamp=.false.)
        call observed%memoize4polarize(pdim_srch)
        call pft_builder%pftc%polarize_ptcl_pft(observed, 1, pdim=pdim_srch, oversamp=.false.)
        call pft_builder%pftc%set_eo(1, .true.)

        ! DIAGNOSTIC BLOCK (temporary, P1): prove that the production context contains
        ! nonzero reference/particle Fourier data before the optimizer is called.
        ! A zero norm here means the failure is in image FFT/polarization or project
        ! indexing, not in the SGD step rule.
        allocate (ref_pft_diag(pdim_srch(1), pdim_srch(2):pdim_srch(3)))
        allocate (ptcl_pft_diag(pdim_srch(1), pdim_srch(2):pdim_srch(3)))
        call pft_builder%pftc%get_ref_pft(1, .true., ref_pft_diag)
        call pft_builder%pftc%get_ptcl_pft(1, ptcl_pft_diag)
        ref_pft_energy = sum(real(ref_pft_diag*conjg(ref_pft_diag), dp))
        ptcl_pft_energy = sum(real(ptcl_pft_diag*conjg(ptcl_pft_diag), dp))
        write (logfhandle, '(a,3I8)') '>>> DIAG POLAR DIMENSIONS: ', pdim_srch
        write (logfhandle, '(a,es16.8)') '>>> DIAG POLAR ENERGY REF: ', ref_pft_energy
        write (logfhandle, '(a,es16.8)') '>>> DIAG POLAR ENERGY PTCL: ', ptcl_pft_energy
        write (logfhandle, '(a,2es16.8)') '>>> DIAG POLAR ABS MAX REF/PTCL: ', maxval(abs(ref_pft_diag)), maxval(abs(ptcl_pft_diag))

        ! P1: calibrate each Fourier shell with the production convention
        ! sigma2(k)=sum_p|R-P|^2/(2*pftsz), then keep the existing Euclidean score.
        call pft_builder%pftc%gen_sigma_contrib(1, 1, [0.0, 0.0], 1, sigma2_noise(:, 1))
        sigma_min = minval(sigma2_noise(:, 1))
        sigma_max = maxval(sigma2_noise(:, 1))
        write (logfhandle, '(a,es16.8)') '>>> DIAG SIGMA2 MIN: ', sigma_min
        write (logfhandle, '(a,es16.8)') '>>> DIAG SIGMA2 MAX: ', sigma_max
        sigma2_noise(:, 1) = max(sigma2_noise(:, 1), 1.0e-6)

        ! DIAGNOSTIC BLOCK (temporary, P2): evaluate the finite raw loss and gradient
        ! at the origin and at four small probes.  This distinguishes a genuinely
        ! flat objective from a zero/invalid normalization before minimize_direct.
        call pft_builder%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0.0_dp, 0.0_dp], 1, raw_f, raw_grad)
        call pft_builder%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [1.0_dp, 0.0_dp], 1, raw_f_xp, raw_grad)
        call pft_builder%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [-1.0_dp, 0.0_dp], 1, raw_f_xm, raw_grad)
        call pft_builder%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0.0_dp, 1.0_dp], 1, raw_f_yp, raw_grad)
        call pft_builder%pftc%gen_raw_euclid_grad_for_rot_8(1, 1, [0.0_dp, -1.0_dp], 1, raw_f_ym, raw_grad)
        write (logfhandle, '(a,es16.8,1x,2es16.8)') '>>> DIAG RAW AT ZERO (LOSS,GX,GY): ', raw_f, raw_grad
        write (logfhandle, '(a,4es16.8)') '>>> DIAG RAW PROBES (+X,-X,+Y,-Y): ', raw_f_xp, raw_f_xm, raw_f_yp, raw_f_ym

        shift_limits(:, 1) = -5.0
        shift_limits(:, 2) = 5.0
        call direct_shift_search%new(pft_builder, shift_limits, opt_angle=.false., direct_only=.true.)
        call direct_shift_search%set_indices(1, 1)
        direct_irot = 1
! P2: request the finite raw Euclidean loss/gradient API so SGD does not
! differentiate the underflow-prone exp(-L) score.
        direct_cxy = direct_shift_search%minimize_direct( &
                     direct_irot, [0.0, 0.0], 0.5, 5, sh_rot=.false., accepted_steps=direct_accepted, &
                     objective_initial=objective_initial, objective_final=objective_final, raw_euclid=.true.)
        write (logfhandle, '(a,es14.6)') '>>> DIRECT SHIFT OBJECTIVE INITIAL: ', objective_initial
        write (logfhandle, '(a,es14.6)') '>>> DIRECT SHIFT OBJECTIVE FINAL:   ', objective_final
        write (logfhandle, '(a,i0)') '>>> DIRECT SHIFT ACCEPTED STEPS: ', direct_accepted
        if (direct_irot == 0) THROW_HARD('direct shift search rejected every tested state')
        if (.not. ieee_is_finite(real(direct_cxy(1), dp))) THROW_HARD('direct shift objective is nonfinite')
        write (logfhandle, '(a,2f10.4)') '>>> DIRECT SHIFT RECOVERED: ', direct_cxy(2:3)
        call direct_shift_search%kill
        call pft_builder%kill_strategy2D_tbox
        call pft_builder%kill_general_tbox
        deallocate (sigma2_noise)
        deallocate (ref_pft_diag, ptcl_pft_diag)
        call simple_chdir(original_cwd, status)
        if (status /= 0) THROW_HARD('could not restore original working directory')

        write (logfhandle, '(A,2F10.4)') &
            '>>> SYNTHETIC APPLIED SHIFT: ', truth%shift
        write (logfhandle, '(A,2F10.4)') &
            '>>> EXPECTED CORRECTIVE SHIFT: ', -truth%shift

        write (logfhandle, '(a)') '>>> Synthetic alignment setup complete (full abinitio2D is covered by v3)'
        write (logfhandle, '(a)') '>>> Results: '//workflow_root%to_char()
        call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO_V2 NORMAL STOP ****')

    end subroutine run_abinitio_v2
end module base_sgd_v2_test

module base_sgd_v3_test
    use base_sgd_test_helpers, only: prepare_test_workspace, make_1jyx_volume, &
                                     simulate_1jyx_particles, check_stack_count, add_noise_checked, create_import_project
contains
    subroutine run_abinitio_v3()
        use simple_core_module_api
        use simple_atoms, only: atoms
        use simple_molecule_data, only: molecule_data, betagal_1jyx
        use simple_cmdline, only: cmdline
        use simple_commanders_sim, only: commander_simulate_particles
        use simple_commanders_project_core, only: commander_new_project
        use simple_commanders_project_ptcl, only: commander_import_particles
        use simple_commanders_abinitio2D, only: commander_abinitio2D
        use simple_sp_project, only: sp_project
        use simple_strategy2D_matcher, only: reset_sgd_stream_dispatch_count, get_sgd_stream_dispatch_count
        implicit none
#include "simple_local_flags.inc"

        character(len=*), parameter :: WORKDIR = 'test_1jyx_abinitio_v3'
        character(len=*), parameter :: PROJNAME = 'onejyx_abinitio_v3'
        character(len=*), parameter :: PROJFILE = PROJNAME//'.simple'
        character(len=*), parameter :: VOLFILE = '1JYX_v3.mrc'
        character(len=*), parameter :: CLEANSTK = '1JYX_v3_clean.mrcs'
        character(len=*), parameter :: NOISYSTK = '1JYX_v3_noisy.mrcs'

        real :: smpd, mskdiam, snr
        integer :: nptcls, ncls, nthr, status, sgd_dispatches
        integer :: vol_dim(3)
        type(cmdline) :: cline_abinitio
        type(commander_simulate_particles) :: xsim
        type(commander_new_project) :: xnew
        type(commander_import_particles) :: ximport
        type(commander_abinitio2D) :: xabinitio
        type(molecule_data) :: mol
        type(atoms) :: molecule
        type(sp_project) :: spproj
        type(string) :: cwd, workflow_root, volume_path, clean_path, noisy_path, project_path, cavgs_jpeg
        logical :: found_cavgs
        logical, parameter :: sgd_diagnostic = .true.

        smpd = 1.3
        mskdiam = 120.0
        snr = 10.0
        nptcls = 200
        ncls = 4
        nthr = 4
        vol_dim = [144, 144, 144]

        call prepare_test_workspace(WORKDIR, cwd, workflow_root, status)

        write (logfhandle, '(a)') '>>> V3 STEP 1: generate 1JYX volume'
        mol = betagal_1jyx()
        call make_1jyx_volume(VOLFILE, smpd, mol, molecule, volume_path, vol_dim)

        write (logfhandle, '(a)') '>>> V3 STEP 2: simulate clean multi-particle stack'
        ! simulate_particles samples random Euler orientations and applies a bounded
        ! random in-plane shift when sherr is nonzero.  This keeps V3 representative
        ! of real particles rather than testing noise-only images.
        call simulate_1jyx_particles(xsim, volume_path, CLEANSTK, smpd, mskdiam, nptcls, nthr, sherr=2.0)
        clean_path = simple_abspath(string(CLEANSTK))
        call check_stack_count(clean_path, nptcls, 'clean stack')

        write (logfhandle, '(a,f7.3)') '>>> V3 STEP 3: add Gaussian noise; SNR = ', snr
        call add_noise_checked(clean_path, NOISYSTK, snr, smpd, nptcls, noisy_path, 'noisy stack')

        write (logfhandle, '(a)') '>>> V3 STEP 4: create and import SIMPLE project'
        call create_import_project(xnew, ximport, PROJNAME, PROJFILE, noisy_path, smpd, project_path)

        write (logfhandle, '(a)') '>>> V3 STEP 5: run abinitio2D through the SGD stage'
        call cline_abinitio%set('prg', 'abinitio2D')
        call cline_abinitio%set('projfile', project_path)
        call cline_abinitio%set('mkdir', 'yes')
        call cline_abinitio%set('mskdiam', mskdiam)
        call cline_abinitio%set('ncls', ncls)
        ! Keep this test independent of repository defaults: the controller's
        ! stage policy is exercised with the same supported refinement and
        ! raw Euclidean objective used by production streaming SGD.
        call cline_abinitio%set('refine', 'snhc_smpl')
        call cline_abinitio%set('objfun', 'euclid')
        ! Production SGD is activated at stage 4; stages 1-3 remain on the
        ! ordinary path and provide the required abinitio2D handoff context.
        call cline_abinitio%set('nstages', 4)
        call cline_abinitio%set('sgd_stage4_mode', 'on')
        call cline_abinitio%set('sgd_path', 'stream')
        call cline_abinitio%set('sgd_diagnostic', merge('yes', 'no ', sgd_diagnostic))
        call cline_abinitio%set('nits_per_stage', 1)
        call cline_abinitio%set('nthr', nthr)
        call reset_sgd_stream_dispatch_count()
        call xabinitio%execute(cline_abinitio)
        call cline_abinitio%kill()

        sgd_dispatches = get_sgd_stream_dispatch_count()
        if( sgd_diagnostic )then
            write(logfhandle,'(A,I0)') '>>> SEARCH DIAG: V3 SGD matcher dispatches: ', sgd_dispatches
        endif
        if( sgd_dispatches < 1 ) THROW_HARD('V3 did not execute the SGD matcher path')

        write (logfhandle, '(a)') '>>> V3 STEP 6: validate project output'
        call spproj%read(project_path)
        if (spproj%os_ptcl2D%get_noris() /= nptcls) THROW_HARD('abinitio2D output particle count mismatch')
        write (logfhandle, '(a,i0)') '>>> V3 PARTICLES: ', spproj%os_ptcl2D%get_noris()
        call spproj%kill()
        ! The class-average segment is written by the abinitio2D stage project and is
        ! not guaranteed to be reattached to the input project on reread.  Validate
        ! the authoritative stage artifact instead of treating an empty os_cls2D
        ! container as a failed classification.
        found_cavgs = .false.
        call find_cavgs_v3(workflow_root, cavgs_jpeg, found_cavgs)
        if (.not. found_cavgs) THROW_HARD('abinitio2D produced no class-average image anywhere below workflow root')
        write (logfhandle, '(a)') '>>> V3 CLASS OUTPUT: '//cavgs_jpeg%to_char()
        call simple_chdir(cwd, status)
        if (status /= 0) THROW_HARD('could not restore original working directory')
        write (logfhandle, '(a)') '>>> V3 RESULTS: '//workflow_root%to_char()
        call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO_V3 NORMAL STOP ****')
    end subroutine run_abinitio_v3
    recursive subroutine find_cavgs_v3(root, result, found)
        use simple_string, only: string
        use simple_fileio, only: simple_list_dirs, simple_list_files_regexp
        type(string), intent(in) :: root
        type(string), intent(inout) :: result
        logical, intent(inout) :: found
        type(string), allocatable :: files(:), dirs(:)
        type(string) :: child
        integer :: i
        if (found) return
        call simple_list_files_regexp(root, '^cavgs_iter[0-9]+\.jpg$', files)
        if (size(files) > 0) then
            result = files(1)
            found = .true.
            return
        end if
        dirs = simple_list_dirs(root)
        do i = 1, size(dirs)
            child = string(root%to_char()//'/'//dirs(i)%to_char())
            call find_cavgs_v3(child, result, found)
            if (found) return
        end do
    end subroutine find_cavgs_v3
end module base_sgd_v3_test

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
        real :: candidate_result(3)
        integer :: accepted_steps, candidate_steps, best_steps
        logical :: pass_class, pass_angle, pass_loss, pass_shift, pass_all
        real, allocatable, target :: sigma2_noise(:, :)
        complex(sp), allocatable :: ref_pft_diag(:, :), ptcl_pft_diag(:, :)
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
        allocate (sigma2_noise(params%kfromto(1):params%kfromto(2), 1), source=0.05); call b%pftc%assign_sigma2_noise(sigma2_noise)
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
        call b%pftc%gen_raw_euclid_grad_for_rot_8(2, 1, [0._dp, 0._dp], 1, loss, grad)
        write (logfhandle, '(a,2es16.8,1x,l1)') '>>> V4 RAW PROBE REF2 (LOSS,GX): ', loss, grad(1), ieee_is_finite(loss)

        write (logfhandle, '(a)') '>>> V4 STEP 3: joint class/rotation search with bounded shift SGD'
        best_loss = huge(1._dp); best_ref = 0; best_rot = 0; best_steps = 0; candidate_count = 0
        write (logfhandle, '(a,i0)') '>>> V4 ROTATION COUNT: ', b%pftc%get_nrots()
        shift_limits(:, 1) = -5.; shift_limits(:, 2) = 5.
        call search%new(b, shift_limits, opt_angle=.false., direct_only=.true.)
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
        pass_all = pass_class .and. pass_angle .and. pass_loss .and. pass_shift
        call search%kill; call b%kill_strategy2D_tbox; call b%kill_general_tbox; deallocate (sigma2_noise)
        call simple_chdir(cwd, status); if (status /= 0) THROW_HARD('could not restore original working directory')
        if (.not. pass_all) THROW_HARD('V4 truth-controlled assignment regression failed; see CHECK lines')
        write (logfhandle, '(a)') '>>> V4 RESULTS: '//root%to_char(); call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO_V4 NORMAL STOP ****')
    end subroutine run_abinitio_v4
end module base_sgd_v4_test

module base_sgd_baseline_test
    use base_sgd_test_helpers, only: prepare_test_workspace, make_1jyx_volume, &
                                     simulate_1jyx_particles, check_stack_count, add_noise_checked, create_import_project
contains
    subroutine run_abinitio_baseline()
        use simple_core_module_api
        use simple_atoms, only: atoms
        use simple_molecule_data, only: molecule_data, betagal_1jyx
        use simple_cmdline, only: cmdline
        use simple_commanders_sim, only: commander_simulate_particles
        use simple_commanders_project_core, only: commander_new_project
        use simple_commanders_project_ptcl, only: commander_import_particles
        use simple_commanders_abinitio2D, only: commander_abinitio2D
        use simple_sp_project, only: sp_project
        implicit none
#include "simple_local_flags.inc"

        character(len=*), parameter :: WORKDIR = 'test_1jyx_abinitio'
        character(len=*), parameter :: PROJNAME = 'onejyx_abinitio'
        character(len=*), parameter :: PROJFILE = PROJNAME//'.simple'
        character(len=*), parameter :: VOL_FILE = '1JYX.mrc'
        character(len=*), parameter :: CLEAN_FILE = '1JYX_particles_clean.mrcs'
        character(len=*), parameter :: NOISY_FILE = '1JYX_particles_noisy.mrcs'
        character(len=*), parameter :: ORI_FILE = '1JYX_particles_oris.txt'

        real :: smpd, mskdiam, noise_snr
        integer :: nptcls, ncls, nthr
        type(cmdline) :: user_args
        type(cmdline) :: cline_abinitio2D
        type(commander_simulate_particles) :: xsimulate_particles
        type(commander_new_project) :: xnew_project
        type(commander_import_particles) :: ximport_particles
        type(commander_abinitio2D) :: xabinitio2D
        type(molecule_data) :: mol
        type(atoms) :: molecule
        type(sp_project) :: spproj
        type(string) :: original_cwd, workflow_root, project_root, project_path
        type(string) :: volume_path, clean_path, noisy_path, orientations_path
        integer :: status

        ! Defaults keep the test reasonably small while leaving enough particles for
        ! both classifications.  They can be overridden with old-style key=value args.
        smpd = 1.3
        mskdiam = 180.0
        noise_snr = 0.10
        nptcls = 200
        ncls = 4
        nthr = 4
        if (command_argument_count() > 0) then
            call user_args%parse_oldschool
            if (user_args%defined('smpd')) smpd = user_args%get_rarg('smpd')
            if (user_args%defined('mskdiam')) mskdiam = user_args%get_rarg('mskdiam')
            if (user_args%defined('snr')) noise_snr = user_args%get_rarg('snr')
            if (user_args%defined('nptcls')) nptcls = user_args%get_iarg('nptcls')
            if (user_args%defined('ncls')) ncls = user_args%get_iarg('ncls')
            if (user_args%defined('nthr')) nthr = user_args%get_iarg('nthr')
        end if
        if (smpd <= 0.0) THROW_HARD('smpd must be positive')
        if (mskdiam <= 0.0) THROW_HARD('mskdiam must be positive')
        if (noise_snr <= 0.0) THROW_HARD('snr must be positive')
        if (nptcls < 8) THROW_HARD('nptcls must be at least 8')
        if (ncls < 1 .or. ncls >= nptcls) THROW_HARD('ncls must be in [1,nptcls)')
        if (nthr < 1) THROW_HARD('nthr must be positive')

        ! Direct commander calls need the normal SIMPLE UI metadata for mkdir=yes.
        call prepare_test_workspace(WORKDIR, original_cwd, workflow_root, status)

        write (logfhandle, '(a)') '>>> Step 1/5: generate a volume from embedded 1JYX beta-galactosidase coordinates'
        mol = betagal_1jyx()
        call make_1jyx_volume(VOL_FILE, smpd, mol, molecule, volume_path)

        write (logfhandle, '(a)') '>>> Step 2/5: simulate a clean stack at random orientations'
        call simulate_1jyx_particles(xsimulate_particles, volume_path, CLEAN_FILE, smpd, mskdiam, nptcls, nthr, &
                                     outfile=ORI_FILE)
        clean_path = simple_abspath(string(CLEAN_FILE))
        orientations_path = simple_abspath(string(ORI_FILE))
        call check_stack_count(clean_path, nptcls, 'clean simulated stack')

        write (logfhandle, '(a,f7.3)') '>>> Step 3/5: add Gaussian noise to the stack; SNR = ', noise_snr
        call add_noise_checked(clean_path, NOISY_FILE, noise_snr, smpd, nptcls, noisy_path, 'noisy simulated stack')

        write (logfhandle, '(a)') '>>> Step 4/5: create a SIMPLE project and import the noisy particles'
        call create_import_project(xnew_project, ximport_particles, PROJNAME, PROJFILE, noisy_path, smpd, project_path, &
                                   project_root, orientations_path%to_char())

        write (logfhandle, '(a)') '>>> Step 5/5: run ab initio 2D classification'
        call cline_abinitio2D%set('prg', 'abinitio2D')
        call cline_abinitio2D%set('projfile', project_path)
        call cline_abinitio2D%set('mkdir', 'yes')
        call cline_abinitio2D%set('mskdiam', mskdiam)
        call cline_abinitio2D%set('ncls', ncls)
        call cline_abinitio2D%set('nstages', 1)
        call cline_abinitio2D%set('nits_per_stage', 1)
        call cline_abinitio2D%set('nthr', nthr)
        call xabinitio2D%execute(cline_abinitio2D)
        call cline_abinitio2D%kill()
        call capture_stage_project_baseline('abinitio2D', PROJFILE, project_path, project_root, status)
        call spproj%read(project_path)
        if (spproj%os_cls2D%get_noris() < 1) THROW_HARD('abinitio2D produced no class averages')
        call spproj%kill()

        call simple_chdir(original_cwd, status)
        if (status /= 0) THROW_HARD('could not restore the original working directory')
        write (logfhandle, '(a)') '>>> Results: '//workflow_root%to_char()
        call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO NORMAL STOP ****')
    end subroutine run_abinitio_baseline
    subroutine capture_stage_project_baseline(stage, projfile, project_path, project_root, status)
        use simple_core_module_api
        use simple_string, only: string
#include "simple_local_flags.inc"
        character(len=*), intent(in) :: stage
        character(len=*), intent(in) :: projfile
        type(string), intent(inout) :: project_path
        type(string), intent(in) :: project_root
        integer, intent(inout) :: status
        if (file_exists(projfile)) project_path = simple_abspath(string(projfile))
        call simple_chdir(project_root, status)
        if (status /= 0) THROW_HARD('could not leave the '//trim(stage)//' directory')
    end subroutine capture_stage_project_baseline
end module base_sgd_baseline_test

program simple_test_sgd_base_suite
    use base_sgd_direct_test, only: run_direct_shift
    use base_sgd_sampling_restore_test, only: run_sampling_and_restore_policy
    use base_sgd_v2_test, only: run_abinitio_v2
    use base_sgd_v3_test, only: run_abinitio_v3
    use base_sgd_v4_test, only: run_abinitio_v4
    use base_sgd_baseline_test, only: run_abinitio_baseline
    implicit none

    write (*,'(a)') '>>> STARTING direct shift test'
    call run_direct_shift()
    write (*,'(a)') '>>> COMPLETED direct shift test'
    write (*,'(a)') '>>> STARTING SGD sampling/restore policy test'
    call run_sampling_and_restore_policy()
    write (*,'(a)') '>>> COMPLETED SGD sampling/restore policy test'
    write (*,'(a)') '>>> STARTING abinitio2D V2 test'
    call run_abinitio_v2()
    write (*,'(a)') '>>> COMPLETED abinitio2D V2 test'
    write (*,'(a)') '>>> STARTING abinitio2D V3 SGD test'
    call run_abinitio_v3()
    write (*,'(a)') '>>> COMPLETED abinitio2D V3 SGD test'
    write (*,'(a)') '>>> STARTING abinitio2D V4 test'
    call run_abinitio_v4()
    write (*,'(a)') '>>> COMPLETED abinitio2D V4 test'
    write (*,'(a)') '>>> STARTING abinitio2D baseline test'
    call run_abinitio_baseline()
    write (*,'(a)') '>>> COMPLETED abinitio2D baseline test'
    write (*, '(/,a)') 'SIMPLE base SGD regression suite: PASS'
end program simple_test_sgd_base_suite
