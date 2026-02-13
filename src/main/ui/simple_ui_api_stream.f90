!@descr: "ui_api_stream" UI api (concrete implementation)
module simple_ui_api_stream
use simple_ui_api_modules
implicit none

type(ui_program), target :: abinitio2D_stream
type(ui_program), target :: assign_optics
type(ui_program), target :: cluster2D_stream
type(ui_program), target :: gen_pickrefs
type(ui_program), target :: pick_extract
type(ui_program), target :: preproc
type(ui_program), target :: sieve_cavgs

contains

    subroutine new_abinitio2D_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio2D_stream%new(&
        &'abinitio2D_stream', &                                                  ! name
        &'2D analysis in streaming mode',&                                       ! descr_short
        &'is a distributed workflow that executes 2D analysis'//&                ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                       ! executable
        &.true.,&                                                                ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,cluster 2D,compute")     ! GUI
        ! image input/output
        ! <empty>
        ! parameter input/output
        call abinitio2D_stream%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', gui_submenu="data", gui_advanced=.false.)
        call abinitio2D_stream%add_input(UI_PARM, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_abinitio2D_stream', .false., '', gui_submenu="data")
        ! alternative inputs
        ! <empty>
        ! search controls
        call abinitio2D_stream%add_input(UI_SRCH, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .true., 200., gui_submenu="cluster 2D",&
        &gui_advanced=.false.)
        ! filter controls
        ! <empty>
        ! mask controls
        call abinitio2D_stream%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., gui_submenu="cluster 2D", gui_advanced=.false.)
        ! computer controls
        call abinitio2D_stream%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call abinitio2D_stream%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        call abinitio2D_stream%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('abinitio2D_stream', abinitio2D_stream, prgtab)
    end subroutine new_abinitio2D_stream


    subroutine new_assign_optics( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call assign_optics%new(&
        &'assign_optics', &                                              ! name
        &'Assign optics groups',&                                        ! descr_short
        &'is a program to assign optics groups during streaming',&       ! descr long
        &'simple_stream',&                                               ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call assign_optics%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls

        call assign_optics%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
    end subroutine new_assign_optics


    subroutine new_cluster2D_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab 
        ! PROGRAM SPECIFICATION
        call cluster2D_stream%new(&
        &'cluster2D_stream', &                                                   ! name
        &'2D analysis in streaming mode',&                                       ! descr_short
        &'is a distributed workflow that executes 2D analysis'//&                ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                       ! executable
        &.true.,&                                                                ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,cluster 2D,compute")     ! GUI
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster2D_stream%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', gui_submenu="data", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_PARM, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_cluster2D_stream', .false., '', gui_submenu="data")
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_stream%add_input(UI_SRCH, ncls_start,                               gui_submenu="cluster 2D", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_SRCH, nptcls_per_cls, required_override=.true., gui_submenu="cluster 2D", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_SRCH, 'ncls', 'num', 'Maximum number of 2D clusters',&
        &'Maximum number of 2D class averages for the pooled particles subsets', 'Maximum # 2D clusters', .true., 200., gui_submenu="cluster 2D",&
        &gui_advanced=.false.)
        ! filter controls
        ! <empty>
        ! mask controls
        call cluster2D_stream%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., gui_submenu="cluster 2D", gui_advanced=.false.)
        ! computer controls
        call cluster2D_stream%add_input(UI_COMP, nchunks,                                gui_submenu="compute", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_COMP, nparts_chunk, required_override=.true., gui_submenu="compute", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_COMP, nparts_pool,  required_override=.true., gui_submenu="compute", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_COMP, nthr,                                   gui_submenu="compute", gui_advanced=.false.)
        call cluster2D_stream%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('cluster2D_stream', cluster2D_stream, prgtab)
    end subroutine new_cluster2D_stream


    subroutine new_gen_pickrefs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call gen_pickrefs%new(&
        &'gen_pickrefs', &                                               ! name
        &'Do a mini stream to create the opening 2D for generation of picking references',&  ! descr_short
        &'is a program to do a mini stream to create the opening 2D',&   ! descr long
        &'simple_stream',&                                               ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! parameter input/output
        call gen_pickrefs%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '')
        call gen_pickrefs%add_input(UI_PARM, 'nmics', 'num', 'Number of micrographs to import',&
        &'Number of micrographs to import for opening 2D', 'Number micrographs', .false., 100.)
        call gen_pickrefs%add_input(UI_PARM, 'optics_dir', 'dir', 'Target directory for optics import',&
        &'Directory where assign_optics application is running', 'e.g. optics_assignment', .false., '')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call gen_pickrefs%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
    end subroutine new_gen_pickrefs


    subroutine new_pick_extract( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call pick_extract%new(&
        &'pick_extract', &                                                               ! name
        &'Preprocessing in streaming mode',&                                             ! descr_short
        &'is a distributed workflow that executes picking and extraction'//&             ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                               ! executable
        &.true.,&                                                                        ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,picking,extract,compute")        ! GUI                    
        ! image input/output
        call pick_extract%add_input(UI_IMG, pickrefs, gui_submenu="picking", gui_advanced=.false., gui_exclusive_group="pickrefs")
        call pick_extract%add_input(UI_IMG, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where a previous pick_extract application was run', 'e.g. 2_pick_extract', .false., '', gui_submenu="data")
        ! parameter input/output
        call pick_extract%add_input(UI_PARM, pcontrast,   gui_submenu="picking")
        call pick_extract%add_input(UI_PARM, box_extract, gui_submenu="extract")
        call pick_extract%add_input(UI_PARM, moldiam,     gui_submenu="picking")
        call pick_extract%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the preprocess_stream application is running', 'e.g. 1_preproc', .true., '', gui_submenu="data")
        call pick_extract%add_input(UI_PARM, 'nmoldiams', 'num', 'Number of molecular diameters to investigate', 'Number of molecular diameters tested',&
        &'e.g. 5', .false., 5., gui_submenu="picking")
        call pick_extract%add_input(UI_PARM, moldiam_max, gui_submenu="picking")
        ! alternative inputs
        ! <empty>
        ! search controls
        call pick_extract%add_input(UI_SRCH, pgrp, required_override=.false., gui_submenu="picking", gui_advanced=.false.)
        ! filter controls
        call pick_extract%add_input(UI_FILT, lp_pick,          gui_submenu="picking")
        call pick_extract%add_input(UI_FILT, ctfresthreshold,  gui_submenu="data")
        call pick_extract%add_input(UI_FILT, icefracthreshold, gui_submenu="data")
        call pick_extract%add_input(UI_FILT, astigthreshold,   gui_submenu="data")
        ! mask controls
        ! <empty>
        ! computer controls
        call pick_extract%add_input(UI_COMP, nthr,   gui_submenu="compute", gui_advanced=.false.)
        call pick_extract%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call pick_extract%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('pick_extract', pick_extract, prgtab)
    end subroutine new_pick_extract

    subroutine new_preproc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call preproc%new(&
        &'preproc', &                                                                       ! name
        &'Preprocessing in streaming mode',&                                                ! descr_short
        &'is a distributed workflow that executes motion_correct, ctf_estimate'//&          ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                                  ! executable
        &.true.,&                                                                           ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,motion correction,CTF estimation")  ! GUI                 
        ! image input/output
        call preproc%add_input(UI_IMG, dir_movies, gui_submenu="data", gui_advanced=.false.)
        call preproc%add_input(UI_IMG, gainref,    gui_submenu="data", gui_advanced=.false.)
        call preproc%add_input(UI_IMG, 'dir_prev', 'file', 'Previous run directory',&
            &'Directory where a previous stream application was run', 'e.g. 2_preproc', .false., '', gui_submenu="data")
        call preproc%add_input(UI_IMG, 'dir_meta', 'dir', 'Directory containing per-movie metadata in XML format',&
            &'Directory containing per-movie metadata XML files from EPU', 'e.g. /dataset/metadata', .false., '', gui_submenu="data", gui_advanced=.false.)
        ! parameter input/output
        call preproc%add_input(UI_PARM, total_dose,                      gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, fraction_dose_target,            gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, smpd_downscale,                  gui_submenu="motion correction", gui_advanced=.false.)
        call preproc%add_input(UI_PARM, eer_fraction,                    gui_submenu="motion correction")
        call preproc%add_input(UI_PARM, max_dose,                        gui_submenu="motion correction")
        call preproc%add_input(UI_PARM, kv,    required_override=.true., gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, cs,    required_override=.true., gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, fraca, required_override=.true., gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, smpd,  required_override=.true., gui_submenu="data",              gui_advanced=.false.)
        call preproc%add_input(UI_PARM, flipgain, gui_submenu="motion correction")
        call preproc%add_input(UI_PARM, 'ninipick', 'num', 'Number of micrographs to perform initial picking preprocessing on',&
        & 'Number of micrographs to perform initial picking preprocessing on', 'e.g 500', .false., 0.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        call preproc%add_input(UI_SRCH, dfmin, gui_submenu="CTF estimation")
        call preproc%add_input(UI_SRCH, dfmax, gui_submenu="CTF estimation")
        call preproc%add_input(UI_SRCH, 'tilt_thres', 'num', 'Threshold for hierarchical clustering of beamtilts',&
        & 'Threshold for hierarchical clustering of beamtilts', 'e.g 0.05', .false., 0.05, gui_submenu="optics groups", gui_online=.true.)
        call preproc%add_input(UI_SRCH, 'beamtilt', 'binary', 'Use beamtilts in optics group assignment',&
        & 'Use beamtilt values (if found in EPU filenames) during optics group assignment(yes|no){yes}', 'beamtilt(yes|no){no}', .false., 'no', gui_submenu="optics groups")
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call preproc%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call preproc%add_input(UI_COMP, nthr,   gui_submenu="compute", gui_advanced=.false.)
        call preproc%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('preproc', preproc, prgtab)
    end subroutine new_preproc

    subroutine new_sieve_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sieve_cavgs%new(&
        &'sieve_cavgs', &                                                       ! name
        &'2D analysis in streaming mode',&                                      ! descr_short
        &'is a distributed workflow that executes 2D analysis'//&               ! descr_long
        &' in streaming mode as the microscope collects the data',&
        &'simple_stream',&                                                      ! executable
        &.true.,&                                                               ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "data,cluster 2D,compute")    ! GUI
        ! image input/output
        ! <empty>
        ! parameter input/output
        call sieve_cavgs%add_input(UI_PARM, 'dir_target', 'file', 'Target directory',&
        &'Directory where the pick_extract application is running', 'e.g. 2_pick_extract', .true., '', gui_submenu="data", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_PARM, 'dir_exec', 'file', 'Previous run directory',&
        &'Directory where previous 2D analysis took place', 'e.g. 3_sieve_cavgs', .false., '', gui_submenu="data")
        ! alternative inputs
        ! <empty>
        ! search controls
        call sieve_cavgs%add_input(UI_SRCH, ncls,                                     gui_submenu="cluster 2D", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_SRCH, nptcls_per_cls, required_override=.true., gui_submenu="cluster 2D", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_SRCH, nchunksperset,                                                      gui_advanced=.false.)
        ! filter controls
        ! <empty>
        ! mask controls
        call sieve_cavgs%add_input(UI_MASK, 'mskdiam', 'num', 'Mask diameter', 'Mask diameter (in A) for application of a soft-edged circular mask to &
        &remove background noise', 'mask diameter in A', .false., 0., gui_submenu="cluster 2D", gui_advanced=.false.)
        ! computer controls
        call sieve_cavgs%add_input(UI_COMP, nchunks,                          gui_submenu="compute", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_COMP, nparts, required_override=.true., gui_submenu="compute", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        call sieve_cavgs%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and management in seconds{1740}(29mins)',&
        &'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('sieve_cavgs', sieve_cavgs, prgtab)
    end subroutine new_sieve_cavgs

end module simple_ui_api_stream
