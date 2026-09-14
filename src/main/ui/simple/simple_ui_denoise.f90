!@descr: module defining the user interfaces for denoising programs in the simple_exec suite
module simple_ui_denoise
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('denoise', 'Denoising', 70)
type(ui_program), target :: icm2D
type(ui_program), target :: icm3D
type(ui_program), target :: ppca_denoise
type(ui_program), target :: ppca_denoise_classes
type(ui_program), target :: denoise_project
type(ui_program), target :: map_params_from_den

contains

    subroutine construct_denoise_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_icm2D(prgtab)
        call new_icm3D(prgtab)
        call new_ppca_denoise(prgtab)
        call new_ppca_denoise_classes(prgtab)
        call new_denoise_project(prgtab)
        call new_map_params_from_den(prgtab)
    end subroutine construct_denoise_programs

    subroutine new_icm2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call icm2D%new(&
        &'icm2D',&                                                                  ! name
        &'Denoise paired 2D image stacks with ICM',& ! summary
        &'is a program for 2D ICM denoising of even/odd image stacks',&             ! help
        &'simple_exec',&                                                            ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call icm2D%add_input(UI_IMG, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call icm2D%add_input(UI_IMG, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call icm2D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call icm2D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call icm2D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('icm2D', icm2D, prgtab, UI_CATEGORY)
    end subroutine new_icm2D

    subroutine new_icm3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call icm3D%new(&
        &'icm3D',&                                                                  ! name
        &'Denoise a 3D volume with nonuniform ICM filtering',& ! summary
        &'is a program for 3D nonuniform filtering by Iterated Conditional Modes',& ! help
        &'simple_exec',&                                                            ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call icm3D%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call icm3D%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call icm3D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call icm3D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! call icm3D%set_input('mask_ctrls', 1, mskdiam)
        ! call icm3D%set_input('mask_ctrls', 2, mskfile)
        ! computer controls
        call icm3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('icm3D', icm3D, prgtab, UI_CATEGORY)
    end subroutine new_icm3D

    subroutine new_ppca_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ppca_denoise%new(&
        &'ppca_denoise',&                             ! name
        &'Denoise an image stack with probabilistic PCA',& ! summary
        &'is a program for ppca-based denoising of an image stack',&  ! help
        &'simple_exec',&                              ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ppca_denoise%add_input(UI_IMG, 'stk',  'file', 'Stack to denoise',  'Stack of images to denoise', 'e.g. stk.mrcs', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call ppca_denoise%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call ppca_denoise%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ppca_denoise%add_input(UI_FILT, 'neigs', 'num', 'Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160)', 'Number of eigencomponents (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs', .true., 160.0, &
        &visibility=UI_VIS_STANDARD)
        call ppca_denoise%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, PPCA plus residual kPCA, standard SVD PCA, kernel PCA, or diffusion maps', 'PCA methods','', .false., 'ppca', &
        &choices=ui_choices([character(len=15) :: 'ppca', 'ppca_kpca_resid', 'pca_svd', 'kpca', 'diffusion_maps']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'k_nn', 'num', 'Diffusion graph neighbors (default 5; try 5-30)', 'Local nearest neighbors used for pca_mode=diffusion_maps', '# neighbors', .false., 5.0, &
        &visibility=UI_VIS_ADVANCED)
            call ppca_denoise%add_input(UI_FILT, 'bandwidth_mode', 'multi', &
                'Diffusion-map bandwidth mode', &
                'Kernel bandwidth policy for diffusion maps(median|ferguson){median}','', .false., 'median', &
            &choices=ui_choices([character(len=8) :: 'median', 'ferguson']), &
        &visibility=UI_VIS_ADVANCED)
            call ppca_denoise%add_input(UI_FILT, 'bandwidth_tune', 'num', &
                'Ferguson bandwidth multiplier (default 1)', &
                'Linear multiplier of the Ferguson-optimal kernel bandwidth (1=optimum); only used when bandwidth_mode=ferguson', &
                'tune >= 0', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
            call ppca_denoise%add_input(UI_FILT, 'dm_alpha', 'num', &
                'Diffusion-map density normalization (default 0)', &
                'Coifman-Lafon alpha: 0=graph Laplacian, 0.5=Fokker-Planck, 1=Laplace-Beltrami (divides out sampling density)', &
                '0 <= alpha <= 1', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}','', .false., 'rbf', &
        &choices=ui_choices([character(len=6) :: 'rbf', 'cosine']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}','', .false., 'nystrom', &
        &choices=ui_choices([character(len=7) :: 'exact', 'nystrom']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'ppca_kpca_resid_alpha', 'num', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'Residual hybrid damping (0 => PPCA only; default 0.5)', 'hybrid alpha', .false., 0.5, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call ppca_denoise%add_input(UI_MASK, mskdiam, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call ppca_denoise%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('ppca_denoise', ppca_denoise, prgtab, UI_CATEGORY)
    end subroutine new_ppca_denoise

    subroutine new_ppca_denoise_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call ppca_denoise_classes%new(&
        &'ppca_denoise_classes',&                     ! name
        &'Denoise class averages with probabilistic PCA',& ! summary
        &'is a program for ppca-based denoising of image classes',&  ! help
        &'simple_exec',&                              ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call ppca_denoise_classes%add_input(UI_PARM, 'pre_norm', 'binary', 'Pre-normalize images', 'Statistical normalization(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call ppca_denoise_classes%add_input(UI_FILT, 'neigs', 'num', '# eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs (0 => auto for Nyström kPCA; default 160; try 128, 160)', '# eigenvecs', .false., 160.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'transp_pca', 'binary', 'transpose for pixel-wise learning', 'transpose for pixel-wise learning(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'pca_mode', 'multi', 'PCA methods: PPCA, standard SVD PCA or kernel PCA', 'PCA methods','', .false., 'ppca', &
        &choices=ui_choices([character(len=7) :: 'ppca', 'pca_svd', 'kpca']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_ker', 'multi', 'Kernel PCA kernel', 'Kernel PCA kernel(rbf|cosine){rbf}','', .false., 'rbf', &
        &choices=ui_choices([character(len=6) :: 'rbf', 'cosine']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_backend', 'multi', 'Kernel PCA backend', 'Kernel PCA backend(exact|nystrom){nystrom}','', .false., 'nystrom', &
        &choices=ui_choices([character(len=7) :: 'exact', 'nystrom']), &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_rbf_gamma', 'num', 'RBF gamma (0 => auto)', 'RBF gamma (0 => auto)', 'gamma', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_nystrom_npts', 'num', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', 'Nyström landmark count (0 => auto=max(128,2*neigs), capped at 512; try 256, 512)', '# landmarks', .false., 512.0, &
        &visibility=UI_VIS_ADVANCED)
        call ppca_denoise_classes%add_input(UI_FILT, 'kpca_nystrom_local_nbrs', 'num', 'Nyström max local support neighbors (default 96; try 96, 128)', 'Nyström max local support neighbors (default 96; try 96, 128)', '# max local nbrs', .false., 96.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call ppca_denoise_classes%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('ppca_denoise_classes', ppca_denoise_classes, prgtab, UI_CATEGORY)
    end subroutine new_ppca_denoise_classes


    subroutine new_denoise_project( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call denoise_project%new(&
        &'denoise_project',&
        &'Create paired raw and denoised particle representations',&
        &'is a workflow for creating a dual-representation project from existing 2D clustering by writing registered phase-flipped raw particles and denoised particle samples from diffusion maps',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'neigs', 'num', &
            'Number of eigencomponents (0 => auto scan; default 200)', &
            'Number of eigencomponents used as the scan upper bound before ICM rank selection', &
            '# eigenvecs', .false., real(DIFFMAP_NEIGS_SCAN_DEFAULT), &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'k_nn', 'num', &
            'Diffusion graph neighbors (default 10; try 5-30)', &
            'Local nearest neighbors used for diffusion-map graph construction', &
            '# neighbors', .false., real(DIFFMAP_GRAPH_KNN_DEFAULT), &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'graph', 'multi', &
            'Diffusion graph', 'Diffusion graph(euc|ori){euc}','', .false., 'euc', &
        &choices=ui_choices([character(len=3) :: 'euc', 'ori']), &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'bandwidth_mode', 'multi', &
            'Diffusion-map bandwidth mode', &
            'Kernel bandwidth policy for diffusion maps(median|ferguson){median}','', .false., 'median', &
        &choices=ui_choices([character(len=8) :: 'median', 'ferguson']), &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'bandwidth_tune', 'num', &
            'Ferguson bandwidth multiplier (default 1)', &
            'Linear multiplier of the Ferguson-optimal kernel bandwidth (1=optimum); only used when bandwidth_mode=ferguson', &
            'tune >= 0', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_FILT, 'dm_alpha', 'num', &
            'Diffusion-map density normalization (default 0)', &
            'Coifman-Lafon alpha: 0=graph Laplacian, 0.5=Fokker-Planck, 1=Laplace-Beltrami (divides out sampling density)', &
            '0 <= alpha <= 1', .false., 0.0, &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_SRCH, nspace, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_SRCH, 'nspace_sub', 'num', &
            'SO3 mixture subspace size', 'SO3 mixture subspace size', &
            '# subspace directions', .false., 500.0, &
        &visibility=UI_VIS_ADVANCED)
        call denoise_project%add_input(UI_MASK, mskdiam, required_override=.false., &
            group="mask", visibility=UI_VIS_STANDARD)
        call denoise_project%add_input(UI_COMP, nparts, required_override=.false., &
            group="compute", visibility=UI_VIS_STANDARD)
        call denoise_project%add_input(UI_COMP, nthr, group="compute", visibility=UI_VIS_STANDARD)
        call add_ui_program('denoise_project', denoise_project, prgtab, UI_CATEGORY)
    end subroutine new_denoise_project

    subroutine new_map_params_from_den( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call map_params_from_den%new(&
        &'map_params_from_den',&
        &'Map denoised-project assignments to raw particles',&
        &'is a workflow for transferring assignments obtained on denoise_project transformed particles back to the raw project particle frame',&
        &'simple_exec',&
        &.true., &
        &visibility=UI_VIS_ADVANCED)
        call map_params_from_den%add_input(UI_FILE, projfile_raw, &
        &visibility=UI_VIS_STANDARD)
        call map_params_from_den%add_input(UI_FILE, projfile_den, &
        &visibility=UI_VIS_STANDARD)
        call map_params_from_den%add_input(UI_FILE, projfile, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call add_ui_program('map_params_from_den', map_params_from_den, prgtab, UI_CATEGORY)
    end subroutine new_map_params_from_den

end module simple_ui_denoise
