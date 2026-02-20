!@descr: module defining the user interfaces for 2D clustering and related programs in the simple_exec suite
module simple_ui_cluster2D
use simple_ui_modules
implicit none

type(ui_program), target :: abinitio2D
type(ui_program), target :: cleanup2D
type(ui_program), target :: cluster2D
type(ui_program), target :: cluster2D_subsets
type(ui_program), target :: make_cavgs
type(ui_program), target :: map_cavgs_selection
type(ui_program), target :: sample_classes
type(ui_program), target :: write_classes

contains

    subroutine construct_cluster2D_programs( prgtab ) 
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio2D(prgtab)
        call new_cleanup2D(prgtab)
        call new_cluster2D_subsets(prgtab)
        call new_make_cavgs(prgtab)
        call new_map_cavgs_selection(prgtab)
        call new_sample_classes(prgtab)
        call new_write_classes(prgtab)
    end subroutine construct_cluster2D_programs

    subroutine print_cluster2D_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('CLUSTER2D WORKFLOWS:', C_UNDERLINED)
        write(logfhandle,'(A)') abinitio2D%name%to_char()
        write(logfhandle,'(A)') cleanup2D%name%to_char()
        write(logfhandle,'(A)') cluster2D%name%to_char()
        write(logfhandle,'(A)') make_cavgs%name%to_char()
        write(logfhandle,'(A)') map_cavgs_selection%name%to_char()
        write(logfhandle,'(A)') sample_classes%name%to_char()
        write(logfhandle,'(A)') write_classes%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_cluster2D_programs

    subroutine new_abinitio2D( prgtab ) 
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio2D%new(&
        &'abinitio2D',&                                                                ! name
        &'ab initio 2D analysis from particles',&                                      ! descr_short
        &'is a distributed workflow for generating 2D class averages from particles',& ! descr_long                                                           ! descr_long
        &'simple_exec',&                                                               ! executable
        &.true.,&                                                                      ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "model,filter,mask,compute"  )       ! GUI
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call abinitio2D%add_input(UI_SRCH, ncls, gui_submenu="search", gui_advanced=.false.)
        call abinitio2D%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their &
        &center of gravity and map shifts back to the particles(yes|no){no}', '(yes|no){no}', .false., 'no', gui_submenu="model")
        call abinitio2D%add_input(UI_SRCH, 'autoscale', 'binary', 'Automatic down-scaling', 'Automatic down-scaling of images &
        &for accelerated computation(yes|no){yes}','(yes|no){yes}', .false., 'yes', gui_submenu="model")
        call abinitio2D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', 'Refinement mode(snhc_smpl|prob|prob_smpl){snhc_smpl}',&
        &'(snhc_smpl|prob|prob_smpl){snhc_smpl}', .false., 'snhc_smpl', gui_submenu="search")
        call abinitio2D%add_input(UI_SRCH, cls_init, gui_submenu="search")
        call abinitio2D%add_input(UI_SRCH, autosample, gui_submenu="search")
        call abinitio2D%add_input(UI_SRCH, 'nsample_start', 'num', 'Starting # of particles per class to sample',&
        &'Starting # of particles per class to sample', 'min # particles per class to sample', .false., 0., gui_submenu="search", gui_advanced=.true.)
        call abinitio2D%add_input(UI_SRCH, 'nsample_stop',  'num', 'Maximum # of particles per class to sample',&
        &'Dynamic particle sampling upper bound to sample', 'max # particles per class to sample', .false., 0., gui_submenu="search", gui_advanced=.true.)
        ! filter controls
        call abinitio2D%add_input(UI_FILT, hp, gui_submenu="filter")
        call abinitio2D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call abinitio2D%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass resolution limit for the first stage of ab-initio model generation',&
            &'low-pass limit in Angstroms', .false., 30., gui_submenu="filter")
        call abinitio2D%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the second stage (no e/o cavgs refinement) in Angstroms', .false., 6., gui_submenu="filter")
        call abinitio2D%add_input(UI_FILT, lp, gui_submenu="filter")
        ! mask controls
        call abinitio2D%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        ! computer controls
        call abinitio2D%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call abinitio2D%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('abinitio2D', abinitio2D, prgtab)
    end subroutine new_abinitio2D

    subroutine new_cleanup2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cleanup2D%new(&
        &'cleanup2D',&                                                           ! name
        &'Analysis of class averages with affinity propagation',&                ! descr_short
        &'is a program for analyzing class averages with affinity propagation',& ! descr_long
        &'simple_exec',&                                                         ! executable
        &.true.)                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! alternative inputs
        ! <empty>
        ! search controls
        ! filter controls
        ! mask controls
        ! computer controls
        ! add 2 ui_hash
        call add_ui_program('cleanup2D', cleanup2D, prgtab)
    end subroutine new_cleanup2D

    subroutine new_cluster2D_subsets( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call cluster2D_subsets%new(&
        &'cluster2D_subsets',&                                                                    ! name
        &'Simultaneous 2D alignment and clustering of single-particle images in streaming mode',& ! descr_short
        &'is a distributed workflow implementing cluster2D in streaming mode',&                   ! descr_long
        &'simple_exec',&                                                                          ! executable
        &.true.,&                                                                                 ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "cluster 2D,compute")                           ! GUI           
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_subsets%add_input(UI_SRCH, nptcls_per_cls, descr_placeholder_override='# of particles per cluster{200}', gui_submenu="cluster 2D", gui_advanced=.false.)
        call cluster2D_subsets%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
            &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes', gui_submenu="cluster 2D")
        call cluster2D_subsets%add_input(UI_SRCH, 'maxnptcls', 'num', 'Maximum # of particles clustered', 'Max # of particles clustered{100000}',&
        &'max # of particles{100000}', .false., 100000., gui_submenu="search", gui_advanced=.true.)
        call cluster2D_subsets%add_input(UI_SRCH, 'nmics', 'num', 'Maximum # of micrographs sampled', 'Max # of micrographs sampled{100}',&
        &'max # of micrographs{100}', .false., 100., gui_submenu="search", gui_advanced=.true.)
        ! filter controls
        call cluster2D_subsets%add_input(UI_FILT, hp, gui_submenu="cluster 2D")
        call cluster2D_subsets%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="cluster 2D")
        call cluster2D_subsets%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Low-pass limit that controls the degree of &
        &downsampling in the second phase. Give estimated best final resolution', 'final low-pass limit in Angstroms', .false., 8.,&
        &gui_submenu="filter", gui_advanced=.true.)
        ! mask controls
        call cluster2D_subsets%add_input(UI_MASK, mskdiam, gui_submenu="cluster 2D", gui_advanced=.false.)
        ! computer controls
        call cluster2D_subsets%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        call cluster2D_subsets%add_input(UI_COMP, 'walltime', 'num', 'Walltime', 'Maximum execution time for job scheduling and &
        &management(29mins){1740}', 'in seconds(29mins){1740}', .false., 1740., gui_submenu="compute")
        ! add to ui_hash
        call add_ui_program('cluster2D_subsets', cluster2D_subsets, prgtab)
    end subroutine new_cluster2D_subsets

    subroutine new_make_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call make_cavgs%new(&
        &'make_cavgs', &                           ! name
        &'Make class averages',&                   ! descr_short
        &'is a distributed workflow for generating class averages or initial random references&
        & for cluster2D execution',&               ! descr_long
        &'simple_exec',&                           ! executable
        &.true.)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call make_cavgs%add_input(UI_IMG, 'refs', 'file', 'Output 2D references',&
        &'Output 2D references', 'xxx.mrc file with references', .false., '')
        ! parameter input/output
        call make_cavgs%add_input(UI_PARM, ncls, required_override=.false.)
        call make_cavgs%add_input(UI_PARM, 'mul', 'num', 'Shift multiplication factor',&
        &'Origin shift multiplication factor{1}','1/scale in pixels{1}', .false., 1.)
        call make_cavgs%add_input(UI_PARM, remap_cls)
        call make_cavgs%add_input(UI_PARM, nspace, required_override=.false.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call make_cavgs%add_input(UI_COMP, nparts)
        call make_cavgs%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('make_cavgs', make_cavgs, prgtab)
    end subroutine new_make_cavgs

    subroutine new_map_cavgs_selection( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call map_cavgs_selection%new(&
        &'map_cavgs_selection',&                                         ! name
        &'Map class average selection to particles in project file',&    ! descr_short
        &'is a program for mapping selection based on class averages to the individual particles using correlation matching',& ! descr_long
        &'all',&                                                         ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call  map_cavgs_selection%add_input(UI_IMG, 'stk', 'file', 'Stack of cavgs to select from', 'Stack of cavgs to select from', 'e.g. cavgs_iter0XX.mrc', .false., '')
        call  map_cavgs_selection%add_input(UI_IMG, 'stk2', 'file', 'Stack of selected cavgs', 'Stack of selected cavgs', 'e.g. selected.spi', .true., '')
        ! parameter input/output
        call  map_cavgs_selection%add_input(UI_PARM, prune)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('map_cavgs_selection', map_cavgs_selection, prgtab)
    end subroutine new_map_cavgs_selection

    subroutine new_sample_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call sample_classes%new(&
        &'sample_classes',&                                                                 ! name
        &'Probabilistic sampling of particles based on class statistics',&                  ! descr_short
        &'is a program for probabilistic sampling of particles based on class statistics',& ! descr_long
        &'simple_exec',&                                                                    ! executable
        &.true.)                                                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call sample_classes%add_input(UI_PARM, 'nptcls_per_part', 'num',    'Number of ptcls per part to select when balancing', '# ptcls per part after balancing', '{100000}', .false., 0.0)
        call sample_classes%add_input(UI_PARM, 'greedy_sampling', 'binary', 'Greedy balanced selection', 'Greedy balanced selection(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call sample_classes%add_input(UI_PARM, 'nparts',          'num',    'Number of partitions in balancing', '# balanced parts', '# balanced parts', .false., 1.)
        call sample_classes%add_input(UI_PARM, nsample)
        call sample_classes%add_input(UI_PARM, 'frac_best',       'num',    'Fraction of best particles to sample from', 'Fraction of best particles to sample from(0-1)', '{0.5}', .false., 0.5)
        call sample_classes%add_input(UI_PARM, 'frac_worst',      'num',    'Fraction of worst particles to sample from', 'Fraction of worst particles to sample from(0-1)', '{0.5}', .false., 0.5)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call sample_classes%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('sample_classes', sample_classes, prgtab)
    end subroutine new_sample_classes

    subroutine new_write_classes( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call write_classes%new(&
        &'write_classes',&                                                                                  ! name
        &'Writes the class averages and the individual (rotated and shifted) particles part of the class',& ! descr_short
        &'is a program for the class averages and the individual (rotated and shifted) particles part of the classto to individual stacks',& ! descr_long
        &'simple_exec',&                                                                                    ! executable
        &.true.)                                                                                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('write_classes', write_classes, prgtab)
    end subroutine new_write_classes

end module simple_ui_cluster2D
