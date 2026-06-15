!@descr: module defining the user interfaces for ab initio 3D reconstruction programs in the simple_exec suite
module simple_ui_abinitio3D
use simple_ui_modules
implicit none

type(ui_program), target :: abinitio3D
type(ui_program), target :: abinitio3D_cavgs
type(ui_program), target :: abinitio3D_cavgs_reject
type(ui_program), target :: estimate_lpstages
type(ui_program), target :: noisevol

contains

    subroutine construct_abinitio3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_abinitio3D(prgtab)
        call new_abinitio3D_cavgs(prgtab)
        call new_abinitio3D_cavgs_reject(prgtab)
        call new_estimate_lpstages(prgtab)
        call new_noisevol(prgtab)
    end subroutine construct_abinitio3D_programs

    subroutine print_abinitio3D_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('AB INITIO 3D RECONSTRUCTION:', C_UNDERLINED)
        write(logfhandle,'(A)') abinitio3D%name%to_char()
        write(logfhandle,'(A)') abinitio3D_cavgs%name%to_char()
        write(logfhandle,'(A)') abinitio3D_cavgs_reject%name%to_char()
        write(logfhandle,'(A)') estimate_lpstages%name%to_char()
        write(logfhandle,'(A)') noisevol%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_abinitio3D_programs

    subroutine new_abinitio3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D%new(&
        &'abinitio3D',&                                                                    ! name
        &'3D ab initio model generation from particles',&                                  ! descr_short
        &'is a distributed workflow for generating an ab initio 3D model from particles',& ! descr_long
        &'simple_exec',&                                                                   ! executable
        &.true.,&                                                                          ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "model,filter,mask,compute"  )           ! GUI
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call abinitio3D%add_input(UI_IMG, 'vol1', 'file', 'Starting template volume', 'Starting reference volume &
        & for particle matching', 'input starting volume e.g. vol.mrc', .false., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call abinitio3D%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){no}', '(yes|no){no}', .false., 'no', gui_submenu="model")
        call abinitio3D%add_input(UI_SRCH, pgrp, gui_submenu="model", gui_advanced=.false.)
        call abinitio3D%add_input(UI_SRCH, pgrp_start, gui_submenu="model")
        call abinitio3D%add_input(UI_SRCH, 'cavg_ini', 'binary', '3D initialization on class averages', '3D initialization on class averages(yes|no){no}', '(yes|no){no}', .false., 'no', gui_submenu="model")
        call abinitio3D%add_input(UI_SRCH, 'cavg_ini_ext', 'binary', 'External class-average 3D initialization', &
            &'Use existing ptcl3D orientations and state assignments from a prior abinitio3D_cavgs run; skips the symmetry-search stage(yes|no){no}', &
            &'(yes|no){no}', .false., 'no', gui_submenu="model", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, nsample, gui_submenu="search", gui_advanced=.false.)
        call abinitio3D%add_input(UI_SRCH, 'nsample_start', 'num', 'Starting number of particles to sample',&
            &'Starting number of particles to sample before ramping to nsample by the stochastic sampling stage; &
            &stage 4 for multivol_mode=independent and stage 5 otherwise',&
            &'starting # particles to sample', .false., 0., gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, 'nstages', 'num', 'Last ab initio stage to run',&
            &'Last abinitio3D stage to run; default is 5 for multivol_mode=independent and 8 otherwise; &
            &independent mode writes final volumes at its last stage',&
            &'last stage', .false., 8., gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, nstates, gui_submenu="search", gui_advanced=.false.)
        call abinitio3D%add_input(UI_SRCH, 'state', 'num', 'Continuation state label', &
            &'State label to select from an existing multi-state abinitio3D project and continue as a single-state stage-5 search', &
            &'state label', .false., 1., gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, 'multivol_mode', 'multi', 'Multi-volume ab initio mode', 'Multi-volume ab initio mode(single|independent|docked){single}', '(single|independent|docked){single}', .false., 'single')
        ! filter controls
        call abinitio3D%add_input(UI_FILT, hp, gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms',  .false., 20., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms; default is 6 for multivol_mode=independent &
            &and 8 otherwise',    .false., 8., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'force_lp_range', 'binary', 'Force low-pass range', &
            &'Use lpstart/lpstop directly for abinitio3D low-pass stages instead of class-FRC-derived limits(yes|no){no}', &
            &'(yes|no){no}', .false., 'no', gui_submenu="filter", gui_advanced=.true.)
        call abinitio3D%add_input(UI_FILT, 'filt_mode', 'multi', 'Filtering mode', &
            &'Filtering mode(none|nonuniform|nonuniform_lpset){nonuniform}; nonuniform_lpset promotes the &
            &NU frontier into an explicit merged-reference LP-set matching run', &
            &'(none|nonuniform|nonuniform_lpset){nonuniform}', .false., 'nonuniform', &
            &gui_submenu="filter", gui_advanced=.true.)
        call abinitio3D%add_input(UI_FILT, conical_fsc, gui_submenu="filter", gui_advanced=.true.)
        call abinitio3D%add_input(UI_FILT, 'lpstart_ini3D',  'num', 'Starting low-pass limit ini3D', 'Starting low-pass limit ini3D',&
            &'low-pass limit for the initial stage of ini3D in Angstroms',  .false., 20., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstop_ini3D',  'num', 'Final low-pass limit ini3D', 'Final low-pass limit ini3D',&
            &'low-pass limit for the final stage of ini3D in Angstroms',    .false., 8., gui_submenu="filter")
        ! mask controls
        call abinitio3D%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        call abinitio3D%add_input(UI_MASK, 'automsk', 'multi', 'Perform envelope masking', &
            &'Whether to generate/apply an envelope mask from the staged automasking point(yes|tight|no){no}', &
            &'(yes|tight|no){no}', .false., 'no', gui_submenu="mask", gui_advanced=.false.)
        ! computer controls
        call abinitio3D%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call abinitio3D%add_input(UI_COMP, nthr,                              gui_submenu="compute", gui_advanced=.false.)
        call abinitio3D%add_input(UI_COMP, 'nthr_ini3D', 'num', 'Number of threads for ini3D phase, give 0 if unsure', 'Number of shared-memory OpenMP threads with close affinity per partition. Typically the same as the number of &
        &logical threads in a socket.', '# shared-memory CPU threads', .false., 0., gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('abinitio3D', abinitio3D, prgtab)
    end subroutine new_abinitio3D

    subroutine new_abinitio3D_cavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D_cavgs%new(&
        &'abinitio3D_cavgs',&                                                                   ! name
        &'3D ab initio model generation from class averages',&                                  ! descr_short
        &'is a distributed workflow for generating an ab initio 3D model from class averages',& ! descr_long
        &'simple_exec',&                                                                        ! executable
        &.true., gui_advanced=.false.)                                                          ! requires sp_project                                         
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call abinitio3D_cavgs%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call abinitio3D_cavgs%add_input(UI_SRCH, pgrp)
        call abinitio3D_cavgs%add_input(UI_SRCH, pgrp_start)
        ! filter controls
        call abinitio3D_cavgs%add_input(UI_FILT, hp, gui_submenu="filter")
        call abinitio3D_cavgs%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call abinitio3D_cavgs%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms', .false., 20., gui_submenu="filter")
        call abinitio3D_cavgs%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms', .false., 8., gui_submenu="filter")
        call abinitio3D_cavgs%add_input(UI_FILT, conical_fsc, gui_submenu="filter", gui_advanced=.true.)
        ! mask controls
        call abinitio3D_cavgs%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        ! computer controls
        call abinitio3D_cavgs%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call abinitio3D_cavgs%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('abinitio3D_cavgs', abinitio3D_cavgs, prgtab)
    end subroutine new_abinitio3D_cavgs

    subroutine new_abinitio3D_cavgs_reject( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call abinitio3D_cavgs_reject%new(&
        &'abinitio3D_cavgs_reject',&                                                                  ! name
        &'Consensus rejection of class averages by restarted multi-state ab initio 3D',&           ! descr_short
        &'runs multiple short two- or three-state abinitio3D_cavgs restarts, builds a consensus state label, &
        &rejects class averages outside the quality-best consensus state, and writes a docked consensus volume',& ! descr_long
        &'simple_exec',&                                                                          ! executable
        &.true.,&                                                                                 ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "model,filter,mask,quality,compute")             ! GUI
        ! search controls
        call abinitio3D_cavgs_reject%add_input(UI_SRCH, 'nstates', 'num', 'Number of ab initio states', &
            &'Number of states used by each short abinitio3D_cavgs restart, either 2 or 3{2}', &
            &'# states{2}', .false., 2., gui_submenu="model", gui_advanced=.false.)
        ! quality controls
        call abinitio3D_cavgs_reject%add_input(UI_SRCH, quality_model, gui_submenu="quality", gui_advanced=.false.)
        call abinitio3D_cavgs_reject%add_input(UI_SRCH, prune, gui_submenu="quality", gui_advanced=.true.)
        call abinitio3D_cavgs_reject%add_input(UI_ALT, 'infile', 'file', 'Quality model input', &
            &'Optional class-average quality model file overriding the built-in preset', &
            &'quality model file', .false., '', gui_submenu="quality", gui_advanced=.true.)
        ! filter controls
        call abinitio3D_cavgs_reject%add_input(UI_FILT, hp, gui_submenu="filter")
        call abinitio3D_cavgs_reject%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call abinitio3D_cavgs_reject%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms', .false., 20., gui_submenu="filter")
        call abinitio3D_cavgs_reject%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms', .false., 8., gui_submenu="filter")
        ! mask controls
        call abinitio3D_cavgs_reject%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        ! compute controls
        call abinitio3D_cavgs_reject%add_input(UI_COMP, 'nrestarts', 'num', 'Number of ab initio restarts', &
            &'Number of asynchronous abinitio3D_cavgs restart jobs to run before consensus voting{3}', &
            &'# restart jobs{3}', .false., 3., gui_submenu="compute", gui_advanced=.false.)
        call abinitio3D_cavgs_reject%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('abinitio3D_cavgs_reject', abinitio3D_cavgs_reject, prgtab)
    end subroutine new_abinitio3D_cavgs_reject

    subroutine new_estimate_lpstages( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call estimate_lpstages%new(&
        &'estimate_lpstages',&                                                                                             ! name
        &'Estimation of low-pass limits, shift boundaries, and downscaling parameters for ab initio 3D',&                  ! descr_short
        &'is a program for estimation of low-pass limits, shift boundaries, and downscaling parameters for ab initio 3D',& ! descr_long
        &'simple_exec',&                                                                                                   ! executable
        &.true.)                                                                                                           ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call estimate_lpstages%add_input(UI_PARM, projfile)
        call estimate_lpstages%add_input(UI_PARM, 'nstages', 'num', 'Number of low-pass limit stages', 'Number of low-pass limit stages', '# stages', .true., 8.)
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
        call add_ui_program('estimate_lpstages', estimate_lpstages, prgtab)
    end subroutine new_estimate_lpstages

    subroutine new_noisevol( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call noisevol%new(&
        &'noisevol',&                         ! name
        &'Generate noise volume',&            ! descr_short
        &'is a program for generating noise volume(s)',&
        &'simple_exec',&                      ! executable
        &.false.)                             ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call noisevol%add_input(UI_PARM, smpd)
        call noisevol%add_input(UI_PARM, box)
        call noisevol%add_input(UI_PARM, 'nstates', 'num', 'Number states', 'Number states', '# states', .false., 1.0)
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
        call add_ui_program('noisevol', noisevol, prgtab)
    end subroutine new_noisevol

end module simple_ui_abinitio3D
