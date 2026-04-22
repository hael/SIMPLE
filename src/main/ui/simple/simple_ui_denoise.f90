!@descr: module defining the user interfaces for denoising programs in the simple_exec suite
module simple_ui_denoise
use simple_ui_modules
implicit none

type(ui_program), target :: icm2D
type(ui_program), target :: icm3D
type(ui_program), target :: ppca_denoise
type(ui_program), target :: ppca_denoise_classes
type(ui_program), target :: ppca_class_splitting
type(ui_program), target :: ppca_volvar

contains

    subroutine construct_denoise_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_icm2D(prgtab)
        call new_icm3D(prgtab)
        call new_ppca_denoise(prgtab)
        call new_ppca_denoise_classes(prgtab)
        call new_ppca_class_splitting(prgtab)
        call new_ppca_volvar(prgtab)
    end subroutine construct_denoise_programs

    subroutine print_denoise_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('DENOISING:', C_UNDERLINED)
        write(logfhandle,'(A)') icm2D%name%to_char()
        write(logfhandle,'(A)') icm3D%name%to_char()
        write(logfhandle,'(A)') ppca_denoise%name%to_char()
        write(logfhandle,'(A)') ppca_denoise_classes%name%to_char()
        write(logfhandle,'(A)') ppca_class_splitting%name%to_char()
        write(logfhandle,'(A)') ppca_volvar%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_denoise_programs

    subroutine new_icm2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! add to ui_hash
        call add_ui_program('icm2D', icm2D, prgtab)
    end subroutine new_icm2D

    subroutine new_icm3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! add to ui_hash
        call add_ui_program('icm3D', icm3D, prgtab)
    end subroutine new_icm3D

    subroutine new_ppca_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ppca_denoise%new(&
        &'ppca_denoise',&                             ! name
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
        call ppca_denoise%add_input(UI_FILT, 'neigs', 'num', 'Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160)', 'Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs', .true., 160.0)
        call ppca_denoise%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, mixture PPCA, PPCA plus residual kPCA, standard SVD PCA or kernel PCA', 'PCA methods', '(ppca|mppca|ppca_kpca_resid|pca_svd|kpca){ppca}', .false., 'ppca')
        call ppca_denoise%add_input(UI_FILT, 'mppca_k', 'num', 'mPPCA mixture components (default 4; try 2, 4, 8)', 'mPPCA mixture components (default 4; try 2, 4, 8)', '# mPPCA comps', .false., 4.0)
        call ppca_denoise%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}', '(rbf|cosine){rbf}', .false., 'rbf')
        call ppca_denoise%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}', '(exact|nystrom){nystrom}', .false., 'nystrom')
        call ppca_denoise%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0)
        call ppca_denoise%add_input(UI_FILT, 'ppca_kpca_resid_alpha', 'num', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'hybrid alpha', .false., 0.5)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0)
        ! mask controls
        ! <empty>
        ! computer controls
        call ppca_denoise%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('ppca_denoise', ppca_denoise, prgtab)
    end subroutine new_ppca_denoise

    subroutine new_ppca_denoise_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        call ppca_denoise_classes%add_input(UI_FILT, 'neigs', 'num', '# eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs', .false., 160.0)
        call ppca_denoise_classes%add_input(UI_FILT, 'transp_pca', 'binary', 'transpose for pixel-wise learning', 'transpose for pixel-wise learning(yes|no){no}', '(yes|no){no}', .false., 'no')
        call ppca_denoise_classes%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, local class-informed PPCA mix, mixture PPCA, PPCA plus residual kPCA, standard SVD PCA or kernel PCA', 'PCA methods', '(ppca|ppca_local_mix|mppca|ppca_kpca_resid|pca_svd|kpca){ppca}', .false., 'ppca')
        call ppca_denoise_classes%add_input(UI_FILT, 'mppca_k', 'num', 'mPPCA mixture components (default 4; try 2, 4, 8)', 'mPPCA mixture components (default 4; try 2, 4, 8)', '# mPPCA comps', .false., 4.0)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}', '(rbf|cosine){rbf}', .false., 'rbf')
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}', '(exact|nystrom){nystrom}', .false., 'nystrom')
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0)
        call ppca_denoise_classes%add_input(UI_FILT, 'ppca_kpca_resid_alpha', 'num', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'hybrid alpha', .false., 0.5)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0)
        ! mask controls
        ! <empty>
        ! computer controls
        call ppca_denoise_classes%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('ppca_denoise_classes', ppca_denoise_classes, prgtab)
    end subroutine new_ppca_denoise_classes

    subroutine new_ppca_volvar( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! add to ui_hash
        call add_ui_program('ppca_volvar', ppca_volvar, prgtab)
    end subroutine new_ppca_volvar

    subroutine new_ppca_class_splitting( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call ppca_class_splitting%new(&
        &'ppca_class_splitting',&
        &'Split classes with shared PPCA latent clustering',&
        &'is a program for splitting 2D/3D particle classes into subclasses using shared PPCA latent spaces and affinity propagation',&
        &'all',&
        &.true.)
        call ppca_class_splitting%add_input(UI_PARM, 'class', 'num', 'Optional class index to split', 'Optional 2D class index or 3D projection/class index to split; omit to process all classes', 'e.g. 5', .false., 0.0)
        call ppca_class_splitting%add_input(UI_PARM, 'pre_norm', 'binary', 'Pre-normalize images', 'Statistical normalization(yes|no){no}', '(yes|no){no}', .false., 'no')
        call ppca_class_splitting%add_input(UI_PARM, 'ncls', 'num', 'Fixed number of subclasses (0 => affinity propagation)', 'Fixed number of subclasses (0 => affinity propagation)', '# subclasses', .false., 0.0)
        call ppca_class_splitting%add_input(UI_PARM, 'nsubcls_max', 'num', 'Maximum subclasses per parent class for AP before k-medoids fallback', 'Maximum subclasses per parent class for AP before k-medoids fallback', '# max subclasses', .false., 8.0)
        call ppca_class_splitting%add_input(UI_ALT,  'oritype', 'multi', 'Particle type to split', 'Particle type to split(ptcl2D|ptcl3D){ptcl2D}', '(ptcl2D|ptcl3D){ptcl2D}', .false., 'ptcl2D')
        call ppca_class_splitting%add_input(UI_FILT, 'neigs', 'num', 'Number of PPCA latent dimensions (0 => auto low-rank scan)', 'Number of PPCA latent dimensions (0 => auto low-rank scan)', '# eigenvecs', .false., 5.0)
        call ppca_class_splitting%add_input(UI_MASK, mskdiam, required_override=.false., gui_submenu="mask", gui_advanced=.false.)
        call ppca_class_splitting%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call ppca_class_splitting%add_input(UI_COMP, nthr,   gui_submenu="compute", gui_advanced=.false.)
        call add_ui_program('ppca_class_splitting', ppca_class_splitting, prgtab)
    end subroutine new_ppca_class_splitting

end module simple_ui_denoise
