!@descr: "single_ui_api_map" UI api (concrete implementation)
module single_ui_api_map
use simple_ui_api_modules
implicit none

type(ui_program), target :: conv_atom_denoise
type(ui_program), target :: tsegmaps_core_finder

contains

    subroutine register_single_ui_map(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('conv_atom_denoise', conv_atom_denoise, prgtab)
        call add_ui_program('tsegmaps_core_finder', tsegmaps_core_finder, prgtab)
    end subroutine register_single_ui_map

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_conv_atom_denoise
        ! PROGRAM SPECIFICATION
        call conv_atom_denoise%new(&
        &'conv_atom_denoise', &                                                  ! name
        &'Denoise atomic-resolution nanoparticle map through atom convolution',& ! descr_short
        &'is a program for denoising atomic-resolution nanoparticle maps exactly as in detect_atoms',& ! descr long
        &'single_exec',&                                                         ! executable
        &.false., gui_advanced=.false.)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call conv_atom_denoise%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Nanoparticle volume to analyse', &
        & 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call conv_atom_denoise%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call conv_atom_denoise%add_input(UI_FILT, element)
        ! mask controls
        call conv_atom_denoise%add_input(UI_MASK, mskdiam)
        ! computer controls
        call conv_atom_denoise%add_input(UI_COMP, nthr)
    end subroutine new_conv_atom_denoise


    subroutine new_tsegmaps_core_finder
        ! PROGRAM SPECIFICATION
        call tsegmaps_core_finder%new(&
        &'tsegmaps_core_finder',&                                                         ! name
        &'For doing radial averaging of the core of docked 3D time-segment maps of NPs',& ! descr_short
        &'is a program that analyses docked time-series density maps',&                   ! descr long
        &'single_exec',&                                                                  ! executable
        &.false., gui_advanced=.false.)                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tsegmaps_core_finder%add_input(UI_IMG, 'filetab', 'file', 'Volumes list',&
        &'List of volumes to analyze', 'list input e.g. voltab.txt', .true., '')
        ! parameter input/output
        call tsegmaps_core_finder%add_input(UI_PARM, smpd)
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
    end subroutine new_tsegmaps_core_finder


end module single_ui_api_map
