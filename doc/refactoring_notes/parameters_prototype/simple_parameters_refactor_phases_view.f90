! View file: semantic structure of the full refactor.
! This is where the current large derivation/checking block gets untangled.

module simple_parameters_refactor_phases_view
implicit none

contains

    subroutine init_execution_context(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors the current EXECUTION RELATED block:
        ! - cwd
        ! - executable
        ! - pid
        ! - distributed worker mode
        ! - program pointer lookup
        ! - previous execution directory detection
    end subroutine init_execution_context

    subroutine attach_program_context(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Determine whether the current program requires a project file,
        ! discover a unique project file when possible, and align prg/projfile/projname.
    end subroutine attach_program_context

    subroutine prepare_execution_directory(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors the mkdir / dir_exec / outdir / project-copying logic.
    end subroutine prepare_execution_directory

    subroutine normalize_path_inputs(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Normalize file and directory paths, including parent-relative adjustments
        ! that currently happen when mkdir=yes.
    end subroutine normalize_path_inputs

    subroutine resolve_project_requirements(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Handle the logic around sp_required, detected *.simple files,
        ! and command-line/project alignment.
    end subroutine resolve_project_requirements

    subroutine resolve_project_inputs(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors the project segment/oritype logic and extraction of
        ! box/smpd/nptcls/ctf/phaseplate from project metadata.
    end subroutine resolve_project_inputs

    subroutine resolve_stack_and_reference_inputs(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Handles stk/oritab/refs-driven defaults and reference file expansion.
    end subroutine resolve_stack_and_reference_inputs

    subroutine resolve_volume_inputs(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Handles vol1 / vol_even / vol_odd / vols(*) consistency and expansion.
    end subroutine resolve_volume_inputs

    subroutine derive_partitioning_defaults(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors split_mode, ncunits, numlen, and basic partition defaults.
    end subroutine derive_partitioning_defaults

    subroutine derive_fractional_update_flags(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors update_frac / trail_rec / fillin / frac_best / frac_worst /
        ! greedy sampling flag derivation.
    end subroutine derive_fractional_update_flags

    subroutine derive_threading_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors nthr/OpenMP handling.
    end subroutine derive_threading_settings

    subroutine validate_walltime(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Keep walltime validation isolated from the rest of parallel setup.
    end subroutine validate_walltime

    subroutine derive_movie_scaling(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors downscale / smpd_downscale / scale_movies.
    end subroutine derive_movie_scaling

    subroutine derive_box_and_sampling(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors box_crop / smpd_crop / xdim / ydim / ldim / clip defaults.
    end subroutine derive_box_and_sampling

    subroutine derive_mask_defaults(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors mskdiam / msk / msk_crop / pftsz logic.
    end subroutine derive_mask_defaults

    subroutine derive_frequency_defaults(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors hp / lp / lpstart / lpstop / fny / dstep / kfromto / lplims2D.
    end subroutine derive_frequency_defaults

    subroutine derive_filter_mode_flags(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors l_lpset / l_envfsc / l_icm / l_gauref / wcrit_enum /
        ! l_graphene / l_potts_prior / l_autoscale / l_sigma_glob /
        ! l_noise_reg / l_lam_anneal / l_ml_reg / l_incrreslim / l_bfac.
    end subroutine derive_filter_mode_flags

    subroutine derive_refine_mode_flags(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Mirrors l_prob_align_mode, neigh defaults, trs defaults,
        ! l_prob_inpl, motion correction convention, and cls_init mode checks.
    end subroutine derive_refine_mode_flags

    subroutine derive_misc_image_flags(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Catch-all for smaller image- and mode-related derived flags that are
        ! not worth mixing into the larger phases above.
    end subroutine derive_misc_image_flags

    subroutine validate_enum_values(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Centralize select-case enum validation for:
        ! - oritype
        ! - objfun
        ! - imgkind
        ! - center_type
        ! - sigma_est
        ! - filt_mode
        ! - mcconvention
        ! - cls_init
    end subroutine validate_enum_values

    subroutine validate_file_combinations(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Centralize conflicting CLI/file combinations such as:
        ! - vol1 vs vol_even/vol_odd
        ! - unsupported mskfile
        ! - refs/ncls consistency
    end subroutine validate_file_combinations

    subroutine validate_dimensions(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! box, box_crop, scale, and other dimension sanity checks.
    end subroutine validate_dimensions

    subroutine validate_image_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Image- and filter-related sanity checks that are easier to maintain
        ! when isolated from derivation.
    end subroutine validate_image_settings

    subroutine validate_refine_settings(self, cline)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        ! Refinement-mode-specific validation and final guardrails.
    end subroutine validate_refine_settings

end module simple_parameters_refactor_phases_view
