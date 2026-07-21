!@descr: module defining the user interfaces for denoising programs in the simple_exec suite
module simple_ui_denoise
use simple_ui_modules
implicit none

type(ui_program), target :: icm2D
type(ui_program), target :: icm3D
type(ui_program), target :: ppca_denoise
type(ui_program), target :: ppca_denoise_classes
type(ui_program), target :: cls_split
type(ui_program), target :: denoise_project
type(ui_program), target :: map_params_from_den
type(ui_program), target :: flex_eigenvol
type(ui_program), target :: ppca_volvar

contains

    subroutine construct_denoise_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_icm2D(prgtab)
        call new_icm3D(prgtab)
        call new_ppca_denoise(prgtab)
        call new_ppca_denoise_classes(prgtab)
        call new_cls_split(prgtab)
        call new_denoise_project(prgtab)
        call new_map_params_from_den(prgtab)
        call new_flex_eigenvol(prgtab)
        call new_ppca_volvar(prgtab)
    end subroutine construct_denoise_programs

    subroutine print_denoise_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('DENOISING:', C_UNDERLINED)
        write(logfhandle,'(A)') icm2D%name%to_char()
        write(logfhandle,'(A)') icm3D%name%to_char()
        write(logfhandle,'(A)') ppca_denoise%name%to_char()
        write(logfhandle,'(A)') ppca_denoise_classes%name%to_char()
        write(logfhandle,'(A)') cls_split%name%to_char()
        write(logfhandle,'(A)') denoise_project%name%to_char()
        write(logfhandle,'(A)') map_params_from_den%name%to_char()
        write(logfhandle,'(A)') flex_eigenvol%name%to_char()
        write(logfhandle,'(A)') ppca_volvar%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_denoise_programs

    subroutine new_icm2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call icm2D%new(&
        &'icm2D',&                                                                  ! name
        &'ICM 2D filter',&                                                          ! descr_short
        &'is a program for 2D ICM denoising of even/odd image stacks',&             ! descr_long
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
        call ppca_denoise%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, PPCA plus residual kPCA, standard SVD PCA, kernel PCA, or diffusion maps', 'PCA methods', '(ppca|ppca_kpca_resid|pca_svd|kpca|diffusion_maps){ppca}', .false., 'ppca')
        call ppca_denoise%add_input(UI_FILT, 'k_nn', 'num', 'Diffusion graph neighbors (default 5; try 5-30)', 'Local nearest neighbors used for pca_mode=diffusion_maps', '# neighbors', .false., 5.0)
        call ppca_denoise%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}', '(rbf|cosine){rbf}', .false., 'rbf')
        call ppca_denoise%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}', '(exact|nystrom){nystrom}', .false., 'nystrom')
        call ppca_denoise%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0)
        call ppca_denoise%add_input(UI_FILT, 'ppca_kpca_resid_alpha', 'num', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'hybrid alpha', .false., 0.5)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0)
        call ppca_denoise%add_input(UI_FILT, hp)
        call ppca_denoise%add_input(UI_FILT, lp)
        ! mask controls
        call ppca_denoise%add_input(UI_MASK, mskdiam, required_override=.false.)
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
        call ppca_denoise_classes%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, standard SVD PCA or kernel PCA', 'PCA methods', '(ppca|pca_svd|kpca){ppca}', .false., 'ppca')
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}', '(rbf|cosine){rbf}', .false., 'rbf')
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}', '(exact|nystrom){nystrom}', .false., 'nystrom')
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0)
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

    subroutine new_cls_split( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call cls_split%new(&
        &'cls_split',&
        &'Split classes with latent clustering',&
        &'splits 2D/3D particle classes into subclasses using diffusion-map or kPCA embeddings and k-medoids clustering',&
        &'all',&
        &.true.)
        call cls_split%add_input(UI_PARM, 'class', 'num', &
            'Optional class index to split', &
            'Optional 2D class index or 3D projection/class index to split; omit to process all classes', &
            'e.g. 5', .false., 0.0)
        call cls_split%add_input(UI_PARM, 'ncls', 'num', &
            'Fixed number of subclasses (0 => auto)', &
            'Fixed number of subclasses (0 => auto)', &
            '# subclasses', .false., 0.0)
        call cls_split%add_input(UI_PARM, 'nsubcls_min', 'num', &
            'Minimum subclass trial count in auto mode (default 3)', &
            'Used only when ncls=0: optimization tries every subclass count from nsubcls_min through nsubcls_max', &
            '# min trial subclasses', .false., 3.0)
        call cls_split%add_input(UI_PARM, 'nsubcls_max', 'num', &
            'Maximum subclass trial count in auto mode (default 10)', &
            'Used only when ncls=0: optimization tries every subclass count from nsubcls_min through nsubcls_max', &
            '# max trial subclasses', .false., 10.0)
        call cls_split%add_input(UI_PARM, 'k_nn', 'num', &
            'Diffusion graph neighbors (default 10; try 5-30)', &
            'Local nearest neighbors used only for diffusion-map modes; larger values smooth the graph', &
            '# neighbors', .false., real(DIFFMAP_GRAPH_KNN_DEFAULT))
        call cls_split%add_input(UI_ALT,  'oritype', 'multi', &
            'Particle type to split', 'Particle type to split(ptcl2D|ptcl3D){ptcl2D}', &
            '(ptcl2D|ptcl3D){ptcl2D}', .false., 'ptcl2D')
        call cls_split%add_input(UI_FILT, 'neigs', 'num', &
            'Number of eigencomponents (0 => auto scan; default 200)', &
            'Number of eigencomponents used as the scan upper bound before ICM dimension selection', &
            '# eigenvecs', .false., real(DIFFMAP_NEIGS_SCAN_DEFAULT))
        call cls_split%add_input(UI_FILT, 'pca_mode', 'multi', &
            'Class split embedding method', &
            'Class split embedding method(diffusion_maps|kpca){diffusion_maps}', &
            '(diffusion_maps|kpca){diffusion_maps}', .false., 'diffusion_maps')
        call cls_split%add_input(UI_FILT, 'graph', 'multi', &
            'Class split graph', 'Class split graph(euc|ori){euc}', &
            '(euc|ori){euc}', .false., 'euc')
        call cls_split%add_input(UI_MASK, mskdiam, required_override=.false., gui_submenu="mask", gui_advanced=.false.)
        call cls_split%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call cls_split%add_input(UI_COMP, nthr,   gui_submenu="compute", gui_advanced=.false.)
        call add_ui_program('cls_split', cls_split, prgtab)
    end subroutine new_cls_split

    subroutine new_denoise_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call denoise_project%new(&
        &'denoise_project',&
        &'Create dual denoised project',&
        &'is a workflow for creating a dual-representation project from existing 2D clustering by writing registered phase-flipped raw particles and denoised particle samples from diffusion maps',&
        &'all',&
        &.true.)
        call denoise_project%add_input(UI_FILT, 'neigs', 'num', &
            'Number of eigencomponents (0 => auto scan; default 200)', &
            'Number of eigencomponents used as the scan upper bound before ICM rank selection', &
            '# eigenvecs', .false., real(DIFFMAP_NEIGS_SCAN_DEFAULT))
        call denoise_project%add_input(UI_FILT, 'k_nn', 'num', &
            'Diffusion graph neighbors (default 10; try 5-30)', &
            'Local nearest neighbors used for diffusion-map graph construction', &
            '# neighbors', .false., real(DIFFMAP_GRAPH_KNN_DEFAULT))
        call denoise_project%add_input(UI_FILT, 'graph', 'multi', &
            'Diffusion graph', 'Diffusion graph(euc|ori){euc}', '(euc|ori){euc}', .false., 'euc')
        call denoise_project%add_input(UI_SRCH, nspace, required_override=.false.)
        call denoise_project%add_input(UI_SRCH, 'nspace_sub', 'num', &
            'SO3 mixture subspace size', 'SO3 mixture subspace size', &
            '# subspace directions', .false., 500.0)
        call denoise_project%add_input(UI_MASK, mskdiam, required_override=.false., &
            gui_submenu="mask", gui_advanced=.false.)
        call denoise_project%add_input(UI_COMP, nparts, required_override=.false., &
            gui_submenu="compute", gui_advanced=.false.)
        call denoise_project%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        call add_ui_program('denoise_project', denoise_project, prgtab)
    end subroutine new_denoise_project

    subroutine new_map_params_from_den( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call map_params_from_den%new(&
        &'map_params_from_den',&
        &'Map denoised-project assignments to raw particles',&
        &'is a workflow for transferring assignments obtained on denoise_project transformed particles back to the raw project particle frame',&
        &'all',&
        &.true.)
        call map_params_from_den%add_input(UI_PARM, projfile_raw)
        call map_params_from_den%add_input(UI_PARM, projfile_den)
        call map_params_from_den%add_input(UI_PARM, projfile, required_override=.false.)
        call add_ui_program('map_params_from_den', map_params_from_den, prgtab)
    end subroutine new_map_params_from_den

    subroutine new_flex_eigenvol( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call flex_eigenvol%new(&
        &'flex_eigenvol',&
        &'Diffusion-manifold representative 3D states',&
        &'builds a sparse diffusion map, selects coverage-initialized k-medoids, and reconstructs soft kernel-weighted 3D pre-image states',&
        &'simple_exec',&
        &.true.)
        call flex_eigenvol%add_input(UI_IMG, 'vol1', 'file', &
            'Mean volume', 'Mean volume used for residual backprojection', &
            'input volume e.g. vol1.mrc', .true., '')
        call flex_eigenvol%add_input(UI_IMG, outvol, required_override=.false.)
        call flex_eigenvol%add_input(UI_SRCH, nspace, required_override=.true.)
        call flex_eigenvol%add_input(UI_FILT, 'neigs', 'num', &
            'Maximum number of diffusion modes (maximum 20)', &
            'Upper bound scanned before ICM feature selection', '# modes', .false., 20.0)
        call flex_eigenvol%add_input(UI_FILT, 'k_nn', 'num', &
            'Nearest neighbors (default 10)', &
            'Registered-residual neighbors retained per particle', '# neighbors', .false., 10.0)
        call flex_eigenvol%add_input(UI_FILT, 'nang_nbrs', 'num', &
            'Angular candidate cap (default 100)', &
            'Maximum orientation-gated candidate particles compared per particle', '# candidates', .false., 100.0)
        call flex_eigenvol%add_input(UI_FILT, 'npreimages', 'num', &
            'Representative state volumes (default 8)', &
            'Number of coverage-initialized k-medoids used as nonlinear diffusion-manifold pre-image targets', &
            '# state volumes', .false., 8.0)
        call flex_eigenvol%add_input(UI_FILT, lp, required_override=.false., &
            descr_placeholder_override='Graph and reconstruction low-pass limit in Angstroms{8}', &
            gui_submenu="regularization", gui_advanced=.false.)
        call flex_eigenvol%add_input(UI_ALT, 'oritype', 'str', &
            'Particle orientation segment', 'Particle orientation segment fixed to ptcl3D', &
            'ptcl3D', .false., 'ptcl3D')
        call flex_eigenvol%add_input(UI_MASK, mskdiam, required_override=.false., &
            gui_submenu="mask", gui_advanced=.false.)
        call flex_eigenvol%add_input(UI_COMP, nparts, required_override=.false., &
            gui_submenu="compute", gui_advanced=.false.)
        call flex_eigenvol%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        call add_ui_program('flex_eigenvol', flex_eigenvol, prgtab)
    end subroutine new_flex_eigenvol

end module simple_ui_denoise
