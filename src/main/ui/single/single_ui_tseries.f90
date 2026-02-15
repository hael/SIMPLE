!@descr: module defining the user interfaces for time-series analysis in the single_exec suite
module single_ui_tseries
use simple_ui_modules
implicit none

type(ui_program), target :: track_particles
type(ui_program), target :: tseries_import
type(ui_program), target :: tseries_make_pickavg
type(ui_program), target :: tseries_motion_correct

contains

    subroutine construct_single_tseries_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_track_particles(prgtab)
        call new_tseries_import(prgtab)
        call new_tseries_make_pickavg(prgtab)
        call new_tseries_motion_correct(prgtab)
    end subroutine construct_single_tseries_programs

    subroutine print_single_tseries_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('TSERIES:', C_UNDERLINED)
        write(logfhandle,'(A)') track_particles%name%to_char()
        write(logfhandle,'(A)') tseries_import%name%to_char()
        write(logfhandle,'(A)') tseries_make_pickavg%name%to_char()
        write(logfhandle,'(A)') tseries_motion_correct%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_tseries_programs

    subroutine new_track_particles( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call track_particles%new(&
        &'track_particles',&                                                     ! name
        &'Track particles in time-series',&                                      ! descr_short
        &'is a distributed workflow for particle tracking in time-series data',& ! descr_long
        &'single_exec',&                                                         ! executable
        &.true., gui_advanced=.false.)                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call track_particles%add_input(UI_PARM, 'fbody', 'string', 'Template output tracked series',&
        &'Template output tracked series', 'e.g. tracked_ptcl', .true., '')
        call track_particles%add_input(UI_PARM, 'boxfile', 'file', 'List of particle coordinates',&
        &'.txt file with EMAN particle coordinates', 'e.g. coords.box', .true., '')
        call track_particles%add_input(UI_PARM, neg)
        call track_particles%add_input(UI_PARM, 'fromf', 'num', 'Frame to start tracking from', 'Frame to start tracking from', 'frame index', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call track_particles%add_input(UI_SRCH, 'offset', 'num', 'Shift half-width search bound', 'Shift half-width search bound(in pixels)',&
        'e.g. pixels window halfwidth', .false., 10.)
        call track_particles%add_input(UI_SRCH, 'nframesgrp', 'num', 'Number of contigous frames to average', '# contigous frames to average before tracking{30}', '{30}', .false., 30.)
        ! <empty>
        ! filter controls
        call track_particles%add_input(UI_FILT, lp_track)
        call track_particles%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the particle and centering', 'centering low-pass limit in Angstroms{5}', .false., 5.)
        call track_particles%add_input(UI_FILT, 'filter', 'multi','Alternative filter for particle tracking',&
            &'Alternative filter for particle tracking(no|tv|nlmean){tv}', '(no|tv|nlmean){tv}', .false., 'tv')
        call track_particles%add_input(UI_FILT, hp)
        ! mask controls
        ! <empty>
        ! computer controls
        call track_particles%add_input(UI_COMP, nthr)
        call add_ui_program('track_particles', track_particles, prgtab)
    end subroutine new_track_particles

    subroutine new_tseries_import( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call tseries_import%new(&
        &'tseries_import',&                                 ! name
        &'Imports time-series datasets',&                   ! descr_short
        &'is a workflow for importing time-series data',&   ! descr_long
        &'single_exec',&                                    ! executable
        &.true., gui_advanced=.false.)                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tseries_import%add_input(UI_IMG, 'filetab', 'file', 'List of individual movie frame files', 'List of frame files (*.mrcs) to import', 'e.g. movie_frames.txt', .true., '')
        ! parameter input/output
        call tseries_import%add_input(UI_PARM, smpd)
        call tseries_import%add_input(UI_PARM, kv, required_override=.true.)
        call tseries_import%add_input(UI_PARM, 'cs', 'num', 'Spherical aberration', 'Spherical aberration constant(in mm){0.0}', 'in mm{0.0}', .true., 0.0)
        call tseries_import%add_input(UI_PARM, 'fraca', 'num', 'Amplitude contrast fraction', 'Fraction of amplitude contrast used for fitting CTF{0.4}', 'fraction{0.4}', .true., 0.4)
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
        call add_ui_program('tseries_import', tseries_import, prgtab)
    end subroutine new_tseries_import

    subroutine new_tseries_make_pickavg( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call tseries_make_pickavg%new(&
        &'tseries_make_pickavg',&                                                        ! name
        &'Align & average the first few frames of the time-series',&                     ! descr_short
        &'is a program for aligning & averaging the first few frames of the time-series&
        & to accomplish SNR enhancement for particle identification',&                   ! descr_long
        &'single_exec',&                                                                 ! executable
        &.true., gui_advanced=.false.)                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_make_pickavg%add_input(UI_PARM, 'nframesgrp', 'num', '# contigous frames to average', 'Number of contigous frames to average using correlation-based weights{10}', '{10}', .false., 10.)
        call tseries_make_pickavg%add_input(UI_PARM, 'fromf',      'num', 'Frame to start averaging from', 'Frame to start averaging from', 'frame index', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_make_pickavg%add_input(UI_SRCH, trs_mc)
        call tseries_make_pickavg%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{5}', .false., 5.)
        call tseries_make_pickavg%add_input(UI_SRCH, mcpatch)
        call tseries_make_pickavg%add_input(UI_SRCH, nxpatch)
        call tseries_make_pickavg%add_input(UI_SRCH, nypatch)
        ! filter controls
        call tseries_make_pickavg%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call tseries_make_pickavg%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){3}', 'in Angstroms{3}', .false., 3.)
        call tseries_make_pickavg%add_input(UI_FILT, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_make_pickavg%add_input(UI_COMP, nthr)
        call add_ui_program('tseries_make_pickavg', tseries_make_pickavg, prgtab)
    end subroutine new_tseries_make_pickavg

    subroutine new_tseries_motion_correct( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call tseries_motion_correct%new(&
        &'tseries_motion_correct', &                                                                               ! name
        &'Anisotropic motion correction of time-series of nanoparticles',&                                         ! descr_short
        &'is a distributed workflow for anisotropic motion correction of time-series (movies) of nanoparticles.',& ! descr_long
        &'single_exec',&                                                                                           ! executable
        &.true., gui_advanced=.false.)                                                                             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call tseries_motion_correct%add_input(UI_PARM, 'boxfile', 'file', 'List of particle coordinates',&
        &'.txt file with EMAN-convention particle coordinates', 'e.g. coords.box', .false., '')        ! alternative inputs
        ! <empty>
        ! search controls
        call tseries_motion_correct%add_input(UI_SRCH, trs_mc)
        call tseries_motion_correct%add_input(UI_SRCH, 'nframesgrp', 'num', '# frames in time moving time window', '# frames in time moving time window subjected to correction', '{5}', .false., 5.)
        call tseries_motion_correct%add_input(UI_SRCH, 'bfac', 'num', 'B-factor applied to frames', 'B-factor applied to frames (in Angstroms^2)', 'in Angstroms^2{5}', .false., 5.)
        call tseries_motion_correct%add_input(UI_SRCH, mcpatch)
        call tseries_motion_correct%add_input(UI_SRCH, nxpatch)
        call tseries_motion_correct%add_input(UI_SRCH, nypatch)
        call tseries_motion_correct%add_input(UI_SRCH, algorithm)
        ! filter controls
        call tseries_motion_correct%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Low-pass limit to be applied in the first &
        &iterations of movie alignment (in Angstroms){5}', 'in Angstroms{5}', .false., 5.)
        call tseries_motion_correct%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit to be applied in the last &
        &iterations of movie alignment (in Angstroms){3}', 'in Angstroms{3}', .false., 3.)
        call tseries_motion_correct%add_input(UI_FILT, wcrit)
        ! mask controls
        ! <empty>
        ! computer controls
        call tseries_motion_correct%add_input(UI_COMP, nparts)
        call tseries_motion_correct%add_input(UI_COMP, nthr)
        call add_ui_program('tseries_motion_correct', tseries_motion_correct, prgtab)
    end subroutine new_tseries_motion_correct

end module single_ui_tseries
