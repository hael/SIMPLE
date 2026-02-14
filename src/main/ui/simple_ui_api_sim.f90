!@descr: "sim" UI api (concrete implementation)
module simple_ui_api_sim
use simple_ui_api_modules
implicit none

type(ui_program), target :: pdb2mrc
type(ui_program), target :: simulate_movie
type(ui_program), target :: simulate_noise
type(ui_program), target :: simulate_particles

contains

    subroutine print_sim_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('SIMULATION:', C_UNDERLINED)
        write(logfhandle,'(A)') pdb2mrc%name%to_char()
        write(logfhandle,'(A)') simulate_movie%name%to_char()
        write(logfhandle,'(A)') simulate_noise%name%to_char()
        write(logfhandle,'(A)') simulate_particles%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_sim_programs

    subroutine new_pdb2mrc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pdb2mrc%new(&
        &'pdb2mrc', &                                      ! name
        &'PDB to MRC simulator',&                          ! descr_short
        &'is a program to simulate a 3D density map in MRC format using a PDB format coordinates file',& ! descr long
        &'all',&                                           ! executable
        &.false., gui_advanced=.false.)                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call pdb2mrc%add_input(UI_IMG, 'pdbfile', 'file', 'PDB input coordinates file', 'Input coordinates file in PDB format', 'PDB file e.g. molecule.pdb', .true., 'molecule.pdb')
        ! parameter input/output
        call pdb2mrc%add_input(UI_PARM, smpd,    required_override=.false.)
        call pdb2mrc%add_input(UI_PARM, vol_dim, required_override=.false.)
        call pdb2mrc%add_input(UI_PARM, outvol)
        call pdb2mrc%add_input(UI_PARM, pdbout)
        call pdb2mrc%add_input(UI_PARM, center_pdb)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('pdb2mrc', pdb2mrc, prgtab)
    end subroutine new_pdb2mrc

    subroutine new_simulate_movie( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_movie%new(&
        &'simulate_movie',&                                 ! name
        &'Simulate DDD movie',&                             ! descr_short
        &'is a program for crude simulation of a DDD movie. Input is a set of projection images to place. &
        &Movie frames are then generated related by randomly shifting the base image and applying noise',& ! descr_long
        &'simple_exec',&                                    ! executable
        &.false.)                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_movie%add_input(UI_IMG, stk)
        ! parameter input/output
        call simulate_movie%add_input(UI_PARM, smpd)
        call simulate_movie%add_input(UI_PARM, 'snr', 'num', 'SNR', 'Signal-to-noise ratio of movie frame', 'signal-to-noise ratio(0.)', .false., 0.)
        call simulate_movie%add_input(UI_PARM, kv)
        call simulate_movie%add_input(UI_PARM, cs)
        call simulate_movie%add_input(UI_PARM, fraca)
        call simulate_movie%add_input(UI_PARM, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2.)
        call simulate_movie%add_input(UI_PARM, trs)
        call simulate_movie%add_input(UI_PARM, 'nframes',  'num', 'Number of frames', 'Number of movie frames', '# frames', .false., 0.)
        call simulate_movie%add_input(UI_PARM, 'xdim',  'num', 'x-dimension', 'Number of pixels in x-direction', '# pixels in x', .false., 0.)
        call simulate_movie%add_input(UI_PARM, 'ydim',  'num', 'y-dimension', 'Number of pixels in y-direction', '# pixels in y', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call simulate_movie%add_input(UI_FILT, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        ! mask controls
        ! <empty>
        ! computer controls
        call simulate_movie%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('simulate_movie', simulate_movie, prgtab)
    end subroutine new_simulate_movie

    subroutine new_simulate_noise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_noise%new(&
        &'simulate_noise',&                                ! name
        &'White noise simulation',&                        ! descr_short
        &'is a program for generating pure noise images',& ! descr_long
        &'simple_exec',&                                   ! executable
        &.false.)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call simulate_noise%add_input(UI_PARM, box)
        call simulate_noise%add_input(UI_PARM, nptcls)
        ! alternative inputs
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
        call add_ui_program('simulate_noise', simulate_noise, prgtab)
    end subroutine new_simulate_noise

    subroutine new_simulate_particles( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call simulate_particles%new(&
        &'simulate_particles',&                                           ! name
        &'Simulate single-particle images',&                              ! descr_short
        &'is a program for simulating single-particle cryo-EM images. It is not a very sophisticated simulator, but&
        & it is nevertheless useful for testing purposes. It does not do any multi-slice simulation and it cannot be&
        & used for simulating molecules containing heavy atoms. It does not even accept a PDB file as an input. Input&
        & is a cryo-EM map, which we usually generate from a PDB file using EMANs program pdb2mrc. The volume is&
        & projected using Fourier interpolation, 20% of the total noise is added to the images (pink noise), they are&
        & then Fourier transformed and multiplied with astigmatic CTF and B-factor. Next, the they are inverse FTed&
        & before the remaining 80% of the noise (white noise) is added',& ! descr_long
        &'simple_exec',&                                                  ! executable
        &.false.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call simulate_particles%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Volume to project', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call simulate_particles%add_input(UI_PARM, smpd)
        call simulate_particles%add_input(UI_PARM, nptcls)
        call simulate_particles%add_input(UI_PARM, 'snr', 'num', 'SNR', 'Signal-to-noise ratio of particle images', 'signal-to-noise ratio(0.)', .true., 0.)
        call simulate_particles%add_input(UI_PARM, oritab)
        call simulate_particles%add_input(UI_PARM, outfile)
        call simulate_particles%add_input(UI_PARM, outstk)
        call simulate_particles%add_input(UI_PARM, 'even', 'binary', 'Generate even projections', 'Generate quasi-even projection directions(yes|no){no}', '(yes|no){no}', .false., 'no')
        call simulate_particles%add_input(UI_PARM, sherr)
        call simulate_particles%add_input(UI_PARM, kv)
        call simulate_particles%add_input(UI_PARM, cs)
        call simulate_particles%add_input(UI_PARM, fraca)
        call simulate_particles%add_input(UI_PARM, deftab)
        call simulate_particles%add_input(UI_PARM, 'defocus',  'num', 'Underfocus', 'Underfocus(in microns)', 'in microns', .false., 2.)
        call simulate_particles%add_input(UI_PARM, dferr)
        call simulate_particles%add_input(UI_PARM, 'astigerr', 'num', 'Astigmatism error', 'Uniform astigmatism error(in microns)', 'error in microns', .false., 0.)
        call simulate_particles%add_input(UI_PARM, ctf)
        call simulate_particles%add_input(UI_PARM, 'nframes', 'num', '# of particle frames', '# of lower SNR particle frames', '{1}', .false., 1.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call simulate_particles%add_input(UI_SRCH, pgrp)
        ! filter controls
        call simulate_particles%add_input(UI_FILT, 'bfac', 'num', 'CTF B-factor','B-factor of CTF in Angstroms^2', 'B-factor in Angstroms^2(>0.0){0}', .false., 0.)
        call simulate_particles%add_input(UI_FILT, 'bfacerr', 'num', 'B-factor error', 'Uniform B-factor error(in Angstroms^2)', 'error(in Angstroms^2)', .false., 50.)
        ! mask controls
        call simulate_particles%add_input(UI_MASK, mskdiam)
        ! computer controls
        call simulate_particles%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('simulate_particles', simulate_particles, prgtab)
    end subroutine new_simulate_particles

end module simple_ui_api_sim
