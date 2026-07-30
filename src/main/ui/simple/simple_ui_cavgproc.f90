!@descr: module defining the user interfaces for class average processing programs in the simple_exec suite
module simple_ui_cavgproc
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('cavgproc', 'Class Average Processing', 40)
type(ui_program), target :: cluster_cavgs
type(ui_program), target :: model_cavgs_rejection
type(ui_program), target :: cluster_cavgs_selection
type(ui_program), target :: cluster_stack
type(ui_program), target :: match_cavgs
type(ui_program), target :: match_stacks
type(ui_program), target :: select_clusters

contains

    subroutine construct_cavgproc_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_cluster_cavgs(prgtab)
        call new_model_cavgs_rejection(prgtab)
        call new_cluster_cavgs_selection(prgtab)
        call new_cluster_stack(prgtab)
        call new_match_cavgs(prgtab)
        call new_match_stacks(prgtab)
        call new_select_clusters(prgtab)
    end subroutine construct_cavgproc_programs

    subroutine print_cavgproc_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('CLASS AVERAGE PROCESSING:', C_UNDERLINED)
        write(logfhandle,'(A)') cluster_cavgs%name%to_char()
        write(logfhandle,'(A)') model_cavgs_rejection%name%to_char()
        write(logfhandle,'(A)') cluster_cavgs_selection%name%to_char()
        write(logfhandle,'(A)') cluster_stack%name%to_char()
        write(logfhandle,'(A)') match_cavgs%name%to_char()
        write(logfhandle,'(A)') match_stacks%name%to_char()
        write(logfhandle,'(A)') select_clusters%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_cavgproc_programs

    subroutine new_cluster_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster_cavgs%new(&
        &'cluster_cavgs',&                                                       ! name
        &'Analysis of class averages with affinity propagation',&                ! summary
        &'is a program for analyzing class averages with affinity propagation',& ! help
        &'simple_exec',&                                                         ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster_cavgs%add_input(UI_IMG, stk, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call cluster_cavgs%add_input(UI_PARM, 'ncls', 'num', 'Number of clusters', 'Number of clusters', '# clusters', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call cluster_cavgs%add_input(UI_PARM, prune, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call cluster_cavgs%add_input(UI_SRCH, clust_crit, &
        &visibility=UI_VIS_ADVANCED)
        call cluster_cavgs%add_input(UI_SRCH, skip_rejection, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call cluster_cavgs%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call cluster_cavgs%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call cluster_cavgs%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call cluster_cavgs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('cluster_cavgs', cluster_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_cluster_cavgs

    subroutine new_model_cavgs_rejection( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call model_cavgs_rejection%new(&
        &'model_cavgs_rejection',&                                               ! name
        &'Model-driven rejection of class averages',&                            ! summary
        &'is a program for automatic class-average rejection using normalized quality feature vectors',& ! help
        &'simple_exec',&                                                         ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                                                 ! requires sp_project except quality_mode=learn|evaluate|promote
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call model_cavgs_rejection%add_input(UI_PARM, quality_mode, &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_PARM, quality_model, &
            activation=ui_activation_equals_any('quality_mode', &
            &[character(len=8) :: 'apply', 'analyze', 'evaluate', 'promote']), &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_PARM, quality_context, &
            activation=ui_activation_equals_any('quality_mode', &
            &[character(len=8) :: 'apply', 'analyze', 'learn', 'evaluate', 'promote']), &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_PARM, prune, &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_FILE, 'filetab', 'file', 'Analysis file table', &
        &'File table of canonical cavgs_quality_training.txt files for quality_mode=learn|evaluate', &
        &'e.g. cavgs_quality_training_filetab.txt', .false., '', &
        &activation=ui_activation_equals_any('quality_mode', [character(len=8) :: 'learn', 'evaluate']), &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_FILE, 'infile', 'file', 'Quality model input', &
        &'Optional learned quality model file for apply/analyze/evaluate or promotion-code generation', &
        &'e.g. cavgs_quality_model_chunk_learned.txt', .false., '', &
        &activation=ui_activation_equals_any('quality_mode', &
        &[character(len=8) :: 'apply', 'analyze', 'evaluate', 'promote']), &
        &visibility=UI_VIS_ADVANCED)
        call model_cavgs_rejection%add_input(UI_FILE, 'fname', 'file', 'Quality model output', &
        &'Output quality model file, evaluation report, or promotion-code snippet for quality_mode=learn|evaluate|promote', &
        &'e.g. cavgs_quality_model_chunk_learned.txt, cavgs_quality_evaluate_report.txt, or cavgs_quality_model_builtin_code.txt', &
        &.false., '', activation=ui_activation_equals_any('quality_mode', &
        &[character(len=8) :: 'learn', 'evaluate', 'promote']), &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call model_cavgs_rejection%add_input(UI_MASK, mskdiam, &
        &activation=ui_activation_equals_any('quality_mode', [character(len=8) :: 'apply', 'analyze', 'evaluate']), &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call model_cavgs_rejection%add_input(UI_COMP, nthr, &
        &activation=ui_activation_equals_any('quality_mode', &
        &[character(len=8) :: 'apply', 'analyze', 'learn', 'evaluate']), &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('model_cavgs_rejection', model_cavgs_rejection, prgtab, UI_CATEGORY)
    end subroutine new_model_cavgs_rejection

    subroutine new_cluster_cavgs_selection( prgtab )    
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster_cavgs_selection%new(&
        &'cluster_cavgs_selection',&                                                                  ! name
        &'Analysis of selected class averages with affinity propagation to prepare for match_cavgs',& ! summary
        &'is a program for analyzing selected class averages with affinity propagation',&             ! help
        &'simple_exec',&                                                                              ! executable
        &.true., &
        &visibility=UI_VIS_ADVANCED)                                                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call cluster_cavgs_selection%add_input(UI_SRCH, clust_crit, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call cluster_cavgs_selection%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call cluster_cavgs_selection%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call cluster_cavgs_selection%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call cluster_cavgs_selection%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('cluster_cavgs_selection', cluster_cavgs_selection, prgtab, UI_CATEGORY)
    end subroutine new_cluster_cavgs_selection

    subroutine new_cluster_stack( prgtab )  
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster_stack%new(&
        &'cluster_stack',&                                            ! name
        &'Analysis of class averages with k-medoids',&                ! summary
        &'is a program for analyzing class averages with k-medoids',& ! help
        &'simple_exec',&                                              ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster_stack%add_input(UI_IMG, stk, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call cluster_stack%add_input(UI_PARM, 'ncls', 'num', 'Number of clusters', 'Number of clusters', '# clusters', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call cluster_stack%add_input(UI_SRCH, clust_crit, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        call cluster_stack%add_input(UI_FILT, hp, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call cluster_stack%add_input(UI_FILT, lp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call cluster_stack%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call cluster_stack%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('cluster_stack', cluster_stack, prgtab, UI_CATEGORY)
    end subroutine new_cluster_stack

    subroutine new_match_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call match_cavgs%new(&
        &'match_cavgs',&                                              ! name
        &'Analysis of class averages with k-medoids',&                ! summary
        &'is a program for analyzing class averages with k-medoids',& ! help
        &'simple_exec',&                                              ! executable
        &.true.)                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call match_cavgs%add_input(UI_FILE, projfile, &
        &visibility=UI_VIS_STANDARD)
        call match_cavgs%add_input(UI_FILE, projfile_ref, &
        &visibility=UI_VIS_STANDARD)
        call match_cavgs%add_input(UI_PARM, prune, &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call match_cavgs%add_input(UI_SRCH, clust_crit, &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call match_cavgs%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_DEVELOPER)
        call match_cavgs%add_input(UI_FILT, lp, &
        &visibility=UI_VIS_DEVELOPER)
        ! mask controls
        call match_cavgs%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call match_cavgs%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('match_cavgs', match_cavgs, prgtab, UI_CATEGORY)
    end subroutine new_match_cavgs

     subroutine new_match_stacks( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call match_stacks%new(&
        &'match_stacks',&                                              ! name
        &'Analysis of class averages with k-medoids',&                ! summary
        &'is a program for analyzing class averages with k-medoids',& ! help
        &'simple_exec',&                                              ! executable
        &.false.)                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call match_stacks%add_input(UI_IMG, stk,  required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call match_stacks%add_input(UI_IMG, stk2, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
         ! <empty>
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call match_stacks%add_input(UI_SRCH, 'clust_crit', 'multi', 'Clustering criterion', 'Clustering criterion(sig|sig_clust|cc|res|hybrid){hybrid}','', .false., 'cc', &
        &choices=ui_choices([character(len=9) :: 'sig', 'sig_clust', 'cc', 'res', 'hybrid']), &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        call match_stacks%add_input(UI_FILT, hp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call match_stacks%add_input(UI_FILT, lp, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call match_stacks%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call match_stacks%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add 2 ui  
        call add_ui_program('match_stacks', match_stacks, prgtab, UI_CATEGORY)
    end subroutine new_match_stacks

    subroutine new_select_clusters( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATIONq
        call select_clusters%new(&
        &'select_clusters',&                                    ! name
        &'Select particle clusters from a SIMPLE project',& ! summary
        &'is a program for selecting clusters from a project',& ! help
        &'simple_exec',&                                        ! executable
        &.true.)                                                ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call select_clusters%add_input(UI_PARM, 'select_flag', 'multi', 'flag to use for selection', 'flag to use for selection (cluster|class){cluster}','', .false., 'cluster', &
        &choices=ui_choices([character(len=7) :: 'cluster', 'class']), &
        &visibility=UI_VIS_DEVELOPER)
        call select_clusters%add_input(UI_PARM, prune, &
        &visibility=UI_VIS_DEVELOPER)
        call select_clusters%add_input(UI_PARM, 'clustinds', 'str', 'Comma separated cluster indices', 'Comma separated cluster indices', 'indx1,indx2', .false., '', &
        &visibility=UI_VIS_DEVELOPER)
        call select_clusters%add_input(UI_PARM, 'clustind',  'num', 'Cluster index', 'Cluster index', 'e.g. 5', .false., 0., &
        &visibility=UI_VIS_DEVELOPER)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('select_clusters', select_clusters, prgtab, UI_CATEGORY)
    end subroutine new_select_clusters

end module simple_ui_cavgproc
