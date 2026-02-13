!@descr: "denoise" UI api (concrete implementation)
module simple_ui_api_denoise
use simple_ui_api_modules
implicit none

type(ui_program), target :: icm2D
type(ui_program), target :: icm3D
type(ui_program), target :: ppca_denoise
type(ui_program), target :: ppca_denoise_classes
type(ui_program), target :: ppca_volvar

contains

    subroutine register_simple_ui_denoise(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('icm2D',                icm2D,                prgtab)
        call add_ui_program('icm3D',                icm3D,                prgtab)
        call add_ui_program('ppca_denoise',         ppca_denoise,         prgtab)
        call add_ui_program('ppca_denoise_classes', ppca_denoise_classes, prgtab)
        call add_ui_program('ppca_volvar',          ppca_volvar,          prgtab)
    end subroutine register_simple_ui_denoise

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_icm2D
        ! PROGRAM SPECIFICATION
        call icm2D%new(&
        &'icm2D',&                                                                  ! name
        &'ICM 2D filter',&                                                          ! descr_short
        &'is a program for 2D nonuniform filtering by Iterated Conditional Modes',& ! descr_long
        &'simple_exec',&                                                            ! executable
        &.false.)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call icm2D%add_input(UI_IMG, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '')
        call icm2D%add_input(UI_IMG, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '')
        ! parameter input/output
        call icm2D%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call icm2D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0)
        ! mask controls
        ! <empty>
        ! computer controls
        call icm2D%add_input(UI_COMP, nthr)
    end subroutine new_icm2D


    subroutine new_icm3D
        ! PROGRAM SPECIFICATION
        call icm3D%new(&
        &'icm3D',&                                                                  ! name
        &'ICM 3D filter',&                                                          ! descr_short
        &'is a program for 3D nonuniform filtering by Iterated Conditional Modes',& ! descr_long
        &'simple_exec',&                                                            ! executable
        &.false.)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call icm3D%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call icm3D%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call icm3D%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call icm3D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0)
        ! mask controls
        ! call icm3D%set_input('mask_ctrls', 1, mskdiam)
        ! call icm3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call icm3D%add_input(UI_COMP, nthr)
    end subroutine new_icm3D


    subroutine new_ppca_denoise
        ! PROGRAM SPECIFICATION
        call ppca_denoise%new(&
        &'ppca_based_denoising',&                     ! name
        &'Filter stack/volume',&                      ! descr_short
        &'is a program for ppca-based denoising of an image stack',&  ! descr_long
        &'simple_exec',&                              ! executable
        &.false.)                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ppca_denoise%add_input(UI_IMG, 'stk',  'file', 'Stack to denoise',  'Stack of images to denoise', 'e.g. stk.mrcs', .true., '')
        call ppca_denoise%add_input(UI_IMG, outstk)
        ! parameter input/output
        call ppca_denoise%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ppca_denoise%add_input(UI_FILT, 'neigs', 'num', 'Number of eigencomponents, corresponding to the number of classes in the stack', 'Number of eigencomponents, corresponding to the number of classes in the stack', '# eigenvecs', .true., 100.0)
        call ppca_denoise%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: probabilistic PCA, standard SVD PCA or kernel PCA', 'PCA methods', '(ppca|pca_svd|kpca){kpca}', .false., 'ppca')
        ! mask controls
        ! <empty>
        ! computer controls
        call ppca_denoise%add_input(UI_COMP, nthr)
    end subroutine new_ppca_denoise


    subroutine new_ppca_denoise_classes
        ! PROGRAM SPECIFICATION
        call ppca_denoise_classes%new(&
        &'ppca_denoise_classes',&                     ! name
        &'Filter stack/volume',&                      ! descr_short
        &'is a program for ppca-based denoising of image classes',&  ! descr_long
        &'all',&                                      ! executable
        &.true.)                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ppca_denoise_classes%add_input(UI_PARM, 'pre_norm', 'binary', 'Pre-normalize images', 'Statistical normalization(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ppca_denoise_classes%add_input(UI_FILT, 'neigs', 'num', '# eigenvecs', '# eigenvecs', '# eigenvecs', .false., 0.0)
        call ppca_denoise_classes%add_input(UI_FILT, 'transp_pca', 'binary', 'transpose for pixel-wise learning', 'transpose for pixel-wise learning(yes|no){no}', '(yes|no){no}', .false., 'no')
        call ppca_denoise_classes%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: probabilistic PCA, standard SVD PCA or kernel PCA', 'PCA methods', '(ppca|pca_svd|kpca){kpca}', .false., 'ppca')
        ! mask controls
        ! <empty>
        ! computer controls
        call ppca_denoise_classes%add_input(UI_COMP, nthr)
    end subroutine new_ppca_denoise_classes


    subroutine new_ppca_volvar
        ! PROGRAM SPECIFICATION
        call ppca_volvar%new(&
        &'ppca_volvar',&                                     ! name
        &'Volume variability analysis using ppca',&          ! descr_short
        &'is a program for ppca-based volume variability',&  ! descr_long
        &'simple_exec',&                                     ! executable
        &.false.)                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ppca_volvar%add_input(UI_IMG, 'vol1', 'file', 'Volume', 'Volume for creating 2D central sections', 'input volume e.g. vol.mrc', .true., 'vol1.mrc')
        ! parameter input/output
        call ppca_volvar%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ppca_volvar%add_input(UI_FILT, 'neigs', 'num', '# eigenvecs', '# eigenvecs', '# eigenvecs', .true., 0.0)
        ! mask controls
        call ppca_volvar%add_input(UI_MASK, mskdiam, required_override=.false.)
        ! computer controls
        call ppca_volvar%add_input(UI_COMP, nthr)
    end subroutine new_ppca_volvar


end module simple_ui_api_denoise
