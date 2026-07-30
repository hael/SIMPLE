!@descr: module defining the user interfaces for map-related programs in the single_exec suite
module single_ui_map
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('map', 'Map Analysis', 60)
type(ui_program), target :: conv_atom_denoise
type(ui_program), target :: tsegmaps_core_finder

contains

    subroutine construct_single_map_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_conv_atom_denoise(prgtab)
        call new_tsegmaps_core_finder(prgtab)
    end subroutine construct_single_map_programs

    subroutine print_single_map_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('MAP ANALYSIS:', C_UNDERLINED)
        write(logfhandle,'(A)') conv_atom_denoise%name%to_char()
        write(logfhandle,'(A)') tsegmaps_core_finder%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_map_programs

    subroutine new_conv_atom_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab           
        ! PROGRAM SPECIFICATION
        call conv_atom_denoise%new(&
        &'conv_atom_denoise', &                                                  ! name
        &'Denoise atomic-resolution nanoparticle map through atom convolution',& ! summary
        &'is a program for denoising atomic-resolution nanoparticle maps exactly as in detect_atoms',& ! descr long
        &'single_exec',&                                                         ! executable
        &.false., visibility=UI_VIS_DEVELOPER)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call conv_atom_denoise%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Nanoparticle volume to analyse', &
        & 'input volume e.g. vol.mrc', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call conv_atom_denoise%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call conv_atom_denoise%add_input(UI_FILT, element, &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call conv_atom_denoise%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call conv_atom_denoise%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('conv_atom_denoise', conv_atom_denoise, prgtab, UI_CATEGORY)
    end subroutine new_conv_atom_denoise

    subroutine new_tsegmaps_core_finder( prgtab )
         class(ui_hash), intent(inout) :: prgtab       
        ! PROGRAM SPECIFICATION
        call tsegmaps_core_finder%new(&
        &'tsegmaps_core_finder',&                                                         ! name
        &'For doing radial averaging of the core of docked 3D time-segment maps of NPs',& ! summary
        &'is a program that analyses docked time-series density maps',&                   ! descr long
        &'single_exec',&                                                                  ! executable
        &.false., visibility=UI_VIS_DEVELOPER)                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call tsegmaps_core_finder%add_input(UI_FILE, 'filetab', 'file', 'Volumes list',&
        &'List of volumes to analyze', 'list input e.g. voltab.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call tsegmaps_core_finder%add_input(UI_PARM, smpd, &
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
        call add_ui_program('tsegmaps_core_finder', tsegmaps_core_finder, prgtab, UI_CATEGORY)
    end subroutine new_tsegmaps_core_finder

end module single_ui_map
