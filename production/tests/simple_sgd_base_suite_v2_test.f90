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
        call direct_shift_search%new_direct(pft_builder, shift_limits)
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
