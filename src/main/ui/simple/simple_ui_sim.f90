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
type(ui_program), target :: simulate_scattering

contains

    subroutine construct_sim_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_cif2mrc(prgtab)
        call new_pdb2mrc(prgtab)
        call new_simulate_movie(prgtab)
        call new_simulate_noise(prgtab)
        call new_simulate_particles(prgtab)
        call new_simulate_scattering(prgtab)
    end subroutine construct_sim_programs

subroutine new_cif2mrc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cif2mrc%new(&
        &'cif2mrc', &                                      ! name
        &'Simulate an MRC density map from PDBx/mmCIF atomic coordinates',& ! summary
        &'is a program to simulate a 3D density map in MRC format using a PDBx/mmCIF format coordinates file',& ! descr long
        &'simple_exec',&                                   ! executable
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
        &'simple_exec',&                                   ! executable
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
        &visibility=UI_VIS_ADVANCED, &
        &display_name='Simulate a Movie')                                           ! requires sp_project
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
        &visibility=UI_VIS_ADVANCED, display_name='Simulate Noise')                                          ! requires sp_project
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
        & it is nevertheless useful for testing purposes. With nstates=K it projects each particle from the volume of&
        & its state (vol1 to volK): the states come from the state field of oritab, or are drawn with the population&
        & fractions of state_pops (equal by default), in exact counts, and are written to the output orientations. Without wave_model it does no multi-slice simulation and it cannot be&
        & used for simulating molecules containing heavy atoms. It does not accept a PDB file as an input. Input&
        & is a cryo-EM map, which we usually generate from a PDB file using EMANs program pdb2mrc. The volume is&
        & projected using Fourier interpolation, 20% of the total noise is added to the images (pink noise), they are&
        & then Fourier transformed and multiplied with astigmatic CTF and B-factor. Next, the they are inverse FTed&
        & before the remaining 80% of the noise (white noise) is added. With wave_model the images are instead those of&
        & the wave model: the map, divided by the scale link of total_dose, is the potential in volts; the predicted&
        & electron counts get Poisson shot noise at the dose, then Gaussian detector noise set by snr',& ! help
        &'simple_exec',&                                                  ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED, &
        &display_name='Simulate Particles')                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_particles%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Volume to project; with nstates=K, vol1 to volK, one per state', 'input volume e.g. vol.mrc', .true., '', &
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
        call simulate_particles%add_input(UI_PARM, wave_model, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, ewald, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, 'total_dose', 'num', 'Electron dose', &
        &'Total exposure in electrons per square Angstrom, with wave_model: the scale link of the map and the&
        & shot noise', 'electrons per Angstrom squared', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, nstates, &
        &visibility=UI_VIS_ADVANCED)
        call simulate_particles%add_input(UI_PARM, 'state_pops', 'str', 'State populations', 'Comma-separated population fractions of the states, summing to 1, every state at least one particle; refused beside an oritab that holds states{equal}', 'e.g. 0.5,0.3,0.2', .false., '', &
        &visibility=UI_VIS_ADVANCED)
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

    subroutine new_simulate_scattering( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call simulate_scattering%new(&
        &'simulate_scattering',&
        &'Simulate wave-scattering cryo-EM images',&
        &'generates detector intensities from a 3D electrostatic potential with the wave models of reconstruct3D_wave (projection, first Born, Rytov on the Ewald paraboloid) or multislice, with Kirkland''s sign conventions and the defocus referred to the specimen midplane. A PDB input is converted with SIMPLE pdb2mrc. The exit wave passes the objective lens of SIMPLE''s CTF before dose and detector noise are applied.',&
        &'simple_exec',&
        &.false., visibility=UI_VIS_ADVANCED, display_name='Simulate Scattering Images')
        call simulate_scattering%add_input(UI_FILE, 'pdbfile', 'file', 'PDB coordinates', &
            &'Atomic coordinates converted to an electrostatic potential with SIMPLE pdb2mrc', &
            &'PDB file, mutually exclusive with vol1', .false., '', visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_IMG, 'vol1', 'file', 'Potential volume', &
            &'Sampled 3D electrostatic potential in volts', 'MRC volume, mutually exclusive with pdbfile', &
            &.false., '', visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_requirement('object', 'Scattering object', &
            &'Select exactly one source for the 3D electrostatic potential', &
            &[character(len=7) :: 'pdbfile', 'vol1'], min_selected=1, max_selected=1)
        call simulate_scattering%add_input(UI_PARM, smpd, required_override=.false., visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, vol_dim, required_override=.false., visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_IMG, outvol, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, nptcls, visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_FILE, oritab, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_FILE, outfile, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_IMG, outstk, visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, wave_model, visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, ewald, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, kv, visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, fraca, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, cs, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, 'defocus', 'num', 'Underfocus', &
            &'Mean underfocus in microns, referred to the specimen midplane (the box centre)', 'microns{2}', .false., 2., &
            &visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, dferr, visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, 'astigerr', 'num', 'Astigmatism error', &
            &'Uniform astigmatism error in microns', 'microns{0}', .false., 0., visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_FILT, 'bfac', 'num', 'Detector envelope B-factor', &
            &'Gaussian Fourier-amplitude envelope B-factor', 'Angstroms squared{64}', .false., 64., &
            &visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_PARM, 'total_dose', 'num', 'Electron dose', &
            &'Poisson exposure dose; zero disables shot noise', 'electrons per Angstrom squared{20}', .false., 20., &
            &visibility=UI_VIS_STANDARD)
        call simulate_scattering%add_input(UI_PARM, 'noise_sdev', 'num', 'Detector noise', &
            &'Additive Gaussian detector-noise standard deviation after dose normalization', &
            &'intensity units{0}', .false., 0., visibility=UI_VIS_ADVANCED)
        call simulate_scattering%add_input(UI_COMP, nthr, visibility=UI_VIS_STANDARD)
        call add_ui_program('simulate_scattering', simulate_scattering, prgtab, UI_CATEGORY)
    end subroutine new_simulate_scattering

end module simple_ui_sim
