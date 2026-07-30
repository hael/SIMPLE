!@descr: module defining the user interfaces for simulation programs in the simple_exec suite
module simple_ui_sim
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('sim', 'Simulation', 140)
type(ui_program), target :: cif2mrc
type(ui_program), target :: pdb2mrc
type(ui_program), target :: simulate_movie
type(ui_program), target :: simulate_noise
type(ui_program), target :: simulate_particles

contains

    subroutine construct_sim_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_cif2mrc(prgtab)
        call new_pdb2mrc(prgtab)
        call new_simulate_movie(prgtab)
        call new_simulate_noise(prgtab)
        call new_simulate_particles(prgtab)
    end subroutine construct_sim_programs

subroutine new_cif2mrc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cif2mrc%new(&
        &'cif2mrc', &                                      ! name
        &'Simulate an MRC density map from PDBx/mmCIF atomic coordinates',& ! summary
        &'is a program to simulate a 3D density map in MRC format using a PDBx/mmCIF format coordinates file',& ! descr long
        &'all',&                                           ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Create Density Map from mmCIF') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cif2mrc%add_input(UI_FILE, 'ciffile', 'file', 'PDBx/mmCIF input coordinates file', 'Input coordinates file in PDBx/mmCIF format', 'PDBx/mmCIF file e.g. molecule.cif', .true., 'molecule.cif', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call cif2mrc%add_input(UI_PARM, smpd,    required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('cif2mrc', cif2mrc, prgtab, UI_CATEGORY)
    end subroutine new_cif2mrc

    subroutine new_pdb2mrc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pdb2mrc%new(&
        &'pdb2mrc', &                                      ! name
        &'Simulate an MRC density map from PDB atomic coordinates',& ! summary
        &'is a program to simulate a 3D density map in MRC format using a PDB format coordinates file',& ! descr long
        &'all',&                                           ! executable
        &.false., visibility=UI_VIS_STANDARD, display_name='Create Density Map from PDB') ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call pdb2mrc%add_input(UI_FILE, 'pdbfile', 'file', 'PDB input coordinates file', 'Input coordinates file in PDB format', 'PDB file e.g. molecule.pdb', .true., 'molecule.pdb', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call pdb2mrc%add_input(UI_PARM, smpd,    required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call pdb2mrc%add_input(UI_PARM, vol_dim, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call pdb2mrc%add_input(UI_IMG, outvol, &
        &visibility=UI_VIS_ADVANCED)
        call pdb2mrc%add_input(UI_FILE, pdbout, &
        &visibility=UI_VIS_ADVANCED)
        call pdb2mrc%add_input(UI_PARM, center_pdb, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('pdb2mrc', pdb2mrc, prgtab, UI_CATEGORY)
    end subroutine new_pdb2mrc

    subroutine new_simulate_movie( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_movie%new(&
        &'simulate_movie',&                                 ! name
        &'Simulate a dose-fractionated electron-microscopy movie',& ! summary
        &'is a program for crude simulation of a DDD movie. Input is a set of projection images to place. &
        &Movie frames are then generated related by randomly shifting the base image and applying noise',& ! help
        &'simple_exec',&                                    ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_movie%add_input(UI_IMG, stk, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call simulate_movie%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call simulate_movie%add_input(UI_PARM, 'snr', 'num', 'SNR', 'Signal-to-noise ratio of movie frame', 'signal-to-noise ratio(0.)', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, kv, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, cs, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, fraca, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, ctf_yes, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, trs, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, 'nframes',  'num', 'Number of frames', 'Number of movie frames', '# frames', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, 'xdim',  'num', 'x-dimension', 'Number of pixels in x-direction', '# pixels in x', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_movie%add_input(UI_PARM, 'ydim',  'num', 'y-dimension', 'Number of pixels in y-direction', '# pixels in y', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call simulate_movie%add_input(UI_FILT, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_movie%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('simulate_movie', simulate_movie, prgtab, UI_CATEGORY)
    end subroutine new_simulate_movie

    subroutine new_simulate_noise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_noise%new(&
        &'simulate_noise',&                                ! name
        &'Generate white-noise images or volumes',& ! summary
        &'is a program for generating pure noise images',& ! help
        &'simple_exec',&                                   ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call simulate_noise%add_input(UI_PARM, box, &
        &visibility=UI_VIS_STANDARD)
        call simulate_noise%add_input(UI_PARM, nptcls, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('simulate_noise', simulate_noise, prgtab, UI_CATEGORY)
    end subroutine new_simulate_noise

    subroutine new_simulate_particles( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_particles%new(&
        &'simulate_particles',&                                           ! name
        &'Simulate single-particle images',&                              ! summary
        &'is a program for simulating single-particle cryo-EM images. It is not a very sophisticated simulator, but&
        & it is nevertheless useful for testing purposes. It does not do any multi-slice simulation and it cannot be&
        & used for simulating molecules containing heavy atoms. It does not even accept a PDB file as an input. Input&
        & is a cryo-EM map, which we usually generate from a PDB file using EMANs program pdb2mrc. The volume is&
        & projected using Fourier interpolation, 20% of the total noise is added to the images (pink noise), they are&
        & then Fourier transformed and multiplied with astigmatic CTF and B-factor. Next, the they are inverse FTed&
        & before the remaining 80% of the noise (white noise) is added',& ! help
        &'simple_exec',&                                                  ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_particles%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Volume to project', 'input volume e.g. vol.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call simulate_particles%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call simulate_particles%add_input(UI_PARM, nptcls, &
        &visibility=UI_VIS_STANDARD)
        call simulate_particles%add_input(UI_PARM, 'snr', 'num', 'SNR', 'Signal-to-noise ratio of particle images', 'signal-to-noise ratio(0.)', .true., 0., &
        &visibility=UI_VIS_STANDARD)
        call simulate_particles%add_input(UI_FILE, oritab, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_FILE, outfile, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, 'even', 'binary', 'Generate even projections', 'Generate quasi-even projection directions(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, sherr, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, kv, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, cs, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, fraca, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_FILE, deftab, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, dferr, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, 'astigerr', 'num', 'Astigmatism error', 'Uniform astigmatism error(in microns)', 'error in microns', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, ctf, &
        &visibility=UI_VIS_STANDARD)
        call simulate_particles%add_input(UI_PARM, 'nframes', 'num', '# of particle frames', '# of lower SNR particle frames', '{1}', .false., 1., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call simulate_particles%add_input(UI_SRCH, pgrp, &
        &visibility=UI_VIS_STANDARD)
        ! filter controls
        call simulate_particles%add_input(UI_FILT, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_FILT, 'bfacerr', 'num', 'B-factor error', 'Uniform B-factor error(in Angstroms^2)', 'error(in Angstroms^2)', .false., 50., &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call simulate_particles%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call simulate_particles%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('simulate_particles', simulate_particles, prgtab, UI_CATEGORY)
    end subroutine new_simulate_particles

end module simple_ui_sim
