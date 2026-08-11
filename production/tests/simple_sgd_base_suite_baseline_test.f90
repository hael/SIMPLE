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
        character(len=256) :: first_arg
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
        integer :: status, first_arg_len, first_arg_status

        ! Defaults keep the test reasonably small while leaving enough particles for
        ! both classifications.  They can be overridden with old-style key=value args.
        smpd = 1.3
        mskdiam = 180.0
        noise_snr = 0.10
        nptcls = 200
        ncls = 4
        nthr = 4
        if (command_argument_count() > 0) then
            call get_command_argument(1, first_arg, first_arg_len, first_arg_status)
            if (first_arg_status /= 0 .or. first_arg_len < 5 .or. first_arg(:5) /= 'case=') then
                call user_args%parse_oldschool
                if (user_args%defined('smpd')) smpd = user_args%get_rarg('smpd')
                if (user_args%defined('mskdiam')) mskdiam = user_args%get_rarg('mskdiam')
                if (user_args%defined('snr')) noise_snr = user_args%get_rarg('snr')
                if (user_args%defined('nptcls')) nptcls = user_args%get_iarg('nptcls')
                if (user_args%defined('ncls')) ncls = user_args%get_iarg('ncls')
                if (user_args%defined('nthr')) nthr = user_args%get_iarg('nthr')
            endif
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
