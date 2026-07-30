!@descr: module defining the user interfaces for masking programs in the simple_exec suite
module simple_ui_mask
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('mask', 'Masking', 100)
type(ui_program), target :: auto_spher_mask
type(ui_program), target :: automask2D
type(ui_program), target :: mask

contains

    subroutine construct_mask_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_auto_spher_mask(prgtab)
        call new_automask2D(prgtab)
        call new_mask(prgtab)
    end subroutine construct_mask_programs

subroutine new_auto_spher_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call auto_spher_mask%new(&
        &'auto_spher_mask',&                              ! name
        &'spherical masking with automatic diameter estimation',& ! summary
        &'is a program for automated spherical masking',& ! help
        &'simple_exec',&                                  ! executable
        &.false.)                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call auto_spher_mask%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call auto_spher_mask%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call auto_spher_mask%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .true., 8., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        ! <empty>
        ! computer controls
        call auto_spher_mask%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('auto_spher_mask', auto_spher_mask, prgtab, UI_CATEGORY)
    end subroutine new_auto_spher_mask

    subroutine new_automask2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call automask2D%new(&
        &'automask2D',&                                        ! name
        &'Create an automatic envelope mask for 2D images',& ! summary
        &'is a program for automated envelope masking in 2D',& ! help
        &'simple_exec',&                                       ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask2D%add_input(UI_IMG, stk, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call automask2D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask2D%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms{20 A}', 'low-pass limit in Angstroms{20 A}', .false., 20., &
        &visibility=UI_VIS_ADVANCED)
        call automask2D%add_input(UI_FILT, 'winsz', 'num', 'Window size for median filter',&
        &'Window size for median filter(in pixels)', 'winsz in pixels', .false., 5.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call automask2D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        call automask2D%add_input(UI_MASK, 'ngrow', 'num', '# layers to grow',&
        &'Binary layers grown for molecular envelope in pixels{3}', 'width of binary layers grown in pixels{3}', .false., 3., &
        &visibility=UI_VIS_ADVANCED)
        call automask2D%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6., &
        &visibility=UI_VIS_ADVANCED)
        call automask2D%add_input(UI_MASK, 'positive', 'binary', 'Consider only positive pixels',&
        &'Consider only positive pixels for threshold determination(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call automask2D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('automask2D', automask2D, prgtab, UI_CATEGORY)
    end subroutine new_automask2D

    subroutine new_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',&                                                        ! name
        &'Apply masks to 2D images or 3D volumes',& ! summary
        &'is a program for masking of 2D images and volumes. If you want to mask your images with a spherical mask with a soft &
        & falloff, set mskdiam to the diameter in A',&                   ! help
        &'simple_exec',&                                                 ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mask%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call mask%add_input(UI_FILE, oritab, &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_FILE, outfile, &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_IMG, stk, &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Volume to mask', &
        & 'input volume e.g. vol.mrc', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_requirement('input', 'Input data', 'Supply either a volume or an image stack.', &
            [character(len=4) :: 'stk ', 'vol1'], max_selected=1)
        ! search controls
        call mask%add_input(UI_SRCH, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity(yes|no){yes}','', .false., 'yes', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call mask%add_input(UI_FILT, lp_backgr, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call mask%add_input(UI_MASK, mskdiam, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_MASK, width, &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6., &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_MASK, 'taper_edges', 'binary', 'Taper edges',&
        &'Whether to taper the edges of image/volume(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call mask%add_input(UI_MASK, 'pdbfile', 'file', 'PDB for 3D envelope masking',&
        &'PDB file used to determine the mask', 'e.g. molecule.pdb', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call mask%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('mask', mask, prgtab, UI_CATEGORY)
    end subroutine new_mask

end module simple_ui_mask
