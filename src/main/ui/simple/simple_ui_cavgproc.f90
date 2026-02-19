!@descr: module defining the user interfaces for class average processing programs in the simple_exec suite
module simple_ui_cavgproc
use simple_ui_modules
implicit none

type(ui_program), target :: cluster_cavgs
type(ui_program), target :: cluster_cavgs_selection
type(ui_program), target :: cluster_stack
type(ui_program), target :: match_cavgs
type(ui_program), target :: match_stacks
type(ui_program), target :: select_clusters

contains

    subroutine construct_cavgproc_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_cluster_cavgs(prgtab)
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
        &'Analysis of class averages with affinity propagation',&                ! descr_short
        &'is a program for analyzing class averages with affinity propagation',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &.true.)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster_cavgs%add_input(UI_PARM, 'ncls', 'num', 'Number of clusters', 'Number of clusters', '# clusters', .false., 0.)
        call cluster_cavgs%add_input(UI_PARM, prune)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster_cavgs%add_input(UI_SRCH, clust_crit)
        ! filter controls
        call cluster_cavgs%add_input(UI_FILT, hp)
        call cluster_cavgs%add_input(UI_FILT, lp)
        ! mask controls
        call cluster_cavgs%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cluster_cavgs%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('cluster_cavgs', cluster_cavgs, prgtab)
    end subroutine new_cluster_cavgs

    subroutine new_cluster_cavgs_selection( prgtab )    
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster_cavgs_selection%new(&
        &'cluster_cavgs_selection',&                                                                  ! name
        &'Analysis of selected class averages with affinity propagation to prepare for match_cavgs',& ! descr_short
        &'is a program for analyzing selected class averages with affinity propagation',&             ! descr_long
        &'simple_exec',&                                                                              ! executable
        &.true.)                                                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster_cavgs_selection%add_input(UI_SRCH, clust_crit)
        ! filter controls
        call cluster_cavgs_selection%add_input(UI_FILT, hp)
        call cluster_cavgs_selection%add_input(UI_FILT, lp)
        ! mask controls
        call cluster_cavgs_selection%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cluster_cavgs_selection%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('cluster_cavgs_selection', cluster_cavgs_selection, prgtab)
    end subroutine new_cluster_cavgs_selection

    subroutine new_cluster_stack( prgtab )  
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster_stack%new(&
        &'cluster_stack',&                                            ! name
        &'Analysis of class averages with k-medoids',&                ! descr_short
        &'is a program for analyzing class averages with k-medoids',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.false.)                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cluster_stack%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call cluster_stack%add_input(UI_PARM, 'ncls', 'num', 'Number of clusters', 'Number of clusters', '# clusters', .false., 0.)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster_stack%add_input(UI_SRCH, clust_crit)
        ! filter controls
        call cluster_stack%add_input(UI_FILT, hp, required_override=.false.)
        call cluster_stack%add_input(UI_FILT, lp, required_override=.true.)
        ! mask controls
        call cluster_stack%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cluster_stack%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('cluster_stack', cluster_stack, prgtab)
    end subroutine new_cluster_stack

    subroutine new_match_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call match_cavgs%new(&
        &'match_cavgs',&                                              ! name
        &'Analysis of class averages with k-medoids',&                ! descr_short
        &'is a program for analyzing class averages with k-medoids',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.true.)                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call match_cavgs%add_input(UI_PARM, projfile)
        call match_cavgs%add_input(UI_PARM, projfile_ref)
        call match_cavgs%add_input(UI_PARM, prune)
        ! alternative inputs
        ! <empty>
        ! search controls
        call match_cavgs%add_input(UI_SRCH, clust_crit)
        ! filter controls
        call match_cavgs%add_input(UI_FILT, hp)
        call match_cavgs%add_input(UI_FILT, lp)
        ! mask controls
        call match_cavgs%add_input(UI_MASK, mskdiam)
        ! computer controls
        call match_cavgs%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('match_cavgs', match_cavgs, prgtab)
    end subroutine new_match_cavgs

     subroutine new_match_stacks( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call match_stacks%new(&
        &'match_stacks',&                                              ! name
        &'Analysis of class averages with k-medoids',&                ! descr_short
        &'is a program for analyzing class averages with k-medoids',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.false.)                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call match_stacks%add_input(UI_IMG, stk,  required_override=.true.)
        call match_stacks%add_input(UI_IMG, stk2, required_override=.true.)
        ! parameter input/output
         ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call match_stacks%add_input(UI_SRCH, 'clust_crit', 'multi', 'Clustering criterion', 'Clustering criterion(sig|sig_clust|cc|res|hybrid){hybrid}',&
        &'(sig|sig_clust|cc|res|hybrid){hybrid}', .false., 'cc')
        ! filter controls
        call match_stacks%add_input(UI_FILT, hp, required_override=.true.)
        call match_stacks%add_input(UI_FILT, lp, required_override=.true.)
        ! mask controls
        call match_stacks%add_input(UI_MASK, mskdiam)
        ! computer controls
        call match_stacks%add_input(UI_COMP, nthr)
        ! add 2 ui  
        call add_ui_program('match_stacks', match_stacks, prgtab)
    end subroutine new_match_stacks

    subroutine new_select_clusters( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATIONq
        call select_clusters%new(&
        &'select_clusters',&                                    ! name
        &'Select clusters',&                                    ! descr_short
        &'is a program for selecting clusters from a project',& ! descr_long
        &'simple_exec',&                                        ! executable
        &.true.)                                                ! requires sp_project
        ! TEMPLATE
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call select_clusters%add_input(UI_PARM, 'select_flag', 'multi', 'flag to use for selection', 'flag to use for selection (cluster|class){cluster}', '(cluster|class){cluster}', .false., 'cluster')
        call select_clusters%add_input(UI_PARM, prune)
        ! alternative inputs
        call select_clusters%add_input(UI_ALT, 'clustinds', 'str', 'Comma separated cluster indices', 'Comma separated cluster indices', 'indx1,indx2', .false., '')
        call select_clusters%add_input(UI_ALT, 'clustind',  'num', 'Cluster index', 'Cluster index', 'e.g. 5', .false., 0.)
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('select_clusters', select_clusters, prgtab)
    end subroutine new_select_clusters

end module simple_ui_cavgproc
