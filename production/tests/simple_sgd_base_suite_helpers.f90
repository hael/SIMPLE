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
