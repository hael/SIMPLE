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
