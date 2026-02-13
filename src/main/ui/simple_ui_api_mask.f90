!@descr: "mask" UI api (concrete implementation)
module simple_ui_api_mask
use simple_ui_api_modules
implicit none

type(ui_program), target :: auto_spher_mask
type(ui_program), target :: automask2D
type(ui_program), target :: mask

contains

    subroutine new_auto_spher_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call auto_spher_mask%new(&
        &'auto_spher_mask',&                              ! name
        &'spherical masking with automatic diameter estimation',& ! descr_short
        &'is a program for automated spherical masking',& ! descr_long
        &'simple_exec',&                                  ! executable
        &.false.)                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call auto_spher_mask%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        ! parameter input/output
        call auto_spher_mask%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call auto_spher_mask%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .true., 8.)
        ! mask controls
        ! <empty>
        ! computer controls
        call auto_spher_mask%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('auto_spher_mask', auto_spher_mask, prgtab) 
    end subroutine new_auto_spher_mask


    subroutine new_automask2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call automask2D%new(&
        &'automask2D',&                                        ! name
        &'2D envelope masking',&                               ! descr_short
        &'is a program for automated envelope masking in 2D',& ! descr_long
        &'simple_exec',&                                       ! executable
        &.false.)                                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask2D%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call automask2D%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask2D%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms{20 A}', 'low-pass limit in Angstroms{20 A}', .false., 20.)
        call automask2D%add_input(UI_FILT, 'winsz', 'num', 'Window size for median filter',&
        &'Window size for median filter(in pixels)', 'winsz in pixels', .false., 5.0)
        ! mask controls
        call automask2D%add_input(UI_MASK, mskdiam)
        call automask2D%add_input(UI_MASK, 'ngrow', 'num', '# layers to grow',&
        &'Binary layers grown for molecular envelope in pixels{3}', 'width of binary layers grown in pixels{3}', .false., 3.)
        call automask2D%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6.)
        call automask2D%add_input(UI_MASK, 'positive', 'binary', 'Consider only positive pixels',&
        &'Consider only positive pixels for threshold determination(yes|no){no}', 'only positive(yes|no){no}', .false., 'no')
        ! computer controls
        call automask2D%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('automask2D', automask2D, prgtab)
    end subroutine new_automask2D

    subroutine new_mask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mask%new(&
        &'mask',&                                                        ! name
        &'Mask images/volumes',&                                         ! descr_short
        &'is a program for masking of 2D images and volumes. If you want to mask your images with a spherical mask with a soft &
        & falloff, set mskdiam to the diameter in A',&                   ! descr_long
        &'simple_exec',&                                                 ! executable
        &.false.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mask%add_input(UI_PARM, smpd)
        call mask%add_input(UI_PARM, oritab)
        call mask%add_input(UI_PARM, outfile)
        ! alternative inputs
        call mask%add_input(UI_ALT, stk)
        call mask%add_input(UI_ALT, 'vol1', 'file', 'Volume', 'Volume to mask', &
        & 'input volume e.g. vol.mrc', .false., '')
        ! search controls
        call mask%add_input(UI_SRCH, 'center', 'binary', 'Center input volume', 'Center input volume by its &
        &center of gravity(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        ! filter controls
        call mask%add_input(UI_FILT, lp_backgr)
        ! mask controls
        call mask%add_input(UI_MASK, mskdiam, required_override=.false.)
        call mask%add_input(UI_MASK, mskfile)
        call mask%add_input(UI_MASK, width)
        call mask%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels', '# pixels cosine edge', .false., 6.)
        call mask%add_input(UI_MASK, 'taper_edges', 'binary', 'Taper edges',&
        &'Whether to taper the edges of image/volume(yes|no){no}', '(yes|no){no}', .false., 'no')
        call mask%add_input(UI_MASK, 'pdbfile', 'file', 'PDB for 3D envelope masking',&
        &'PDB file used to determine the mask', 'e.g. molecule.pdb', .false., '')
        ! computer controls
        call mask%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('mask', mask, prgtab)
    end subroutine new_mask

end module simple_ui_api_mask
