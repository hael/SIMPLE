module simple_ui_api_abinitio3D
use simple_ui_api_modules
implicit none

type(ui_program), target :: abinitio3D
type(ui_program), target :: abinitio3D_cavgs
type(ui_program), target :: estimate_lpstages
type(ui_program), target :: multivol_assign
type(ui_program), target :: noisevol

contains

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
        ! <empty>
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
        call abinitio3D%add_input(UI_SRCH, nsample, gui_submenu="search", gui_advanced=.false.)
        call abinitio3D%add_input(UI_SRCH, 'nsample_start', 'num', 'Dynamic particle sampling lower bound', 'Dynamic particle sampling lower bound', 'min # particles to sample', .false., 0., gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, 'nsample_stop',  'num', 'Dynamic particle sampling upper bound', 'Dynamic particle sampling upper bound', 'max # particles to sample', .false., 0., gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, update_frac, gui_submenu="search", gui_advanced=.true.)
        call abinitio3D%add_input(UI_SRCH, nstates, gui_submenu="search", gui_advanced=.false.)
        call abinitio3D%add_input(UI_SRCH, 'multivol_mode', 'multi', 'Multi-volume ab initio mode', 'Multi-volume ab initio mode(single|independent|docked){single}', '(single|independent|docked){single}', .false., 'single')
        ! filter controls
        call abinitio3D%add_input(UI_FILT, hp, gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms',  .false., 20., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms',    .false., 8., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstart_ini3D',  'num', 'Starting low-pass limit ini3D', 'Starting low-pass limit ini3D',&
            &'low-pass limit for the initial stage of ini3D in Angstroms',  .false., 20., gui_submenu="filter")
        call abinitio3D%add_input(UI_FILT, 'lpstop_ini3D',  'num', 'Final low-pass limit ini3D', 'Final low-pass limit ini3D',&
            &'low-pass limit for the final stage of ini3D in Angstroms',    .false., 8., gui_submenu="filter")
        ! mask controls
        call abinitio3D%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
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
        ! mask controls
        call abinitio3D_cavgs%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        ! computer controls
        call abinitio3D_cavgs%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('abinitio3D_cavgs', abinitio3D_cavgs, prgtab)
    end subroutine new_abinitio3D_cavgs


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


    subroutine new_multivol_assign( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call multivol_assign%new(&
        &'multivol_assign',&                                                               ! name
        &'multi-volume assignment and 3D reconstruction from particles',&                  ! descr_short
        &'is a distributed workflow for generating multiple structural state volumes from particles',& ! descr_long                                                         ! descr_long
        &'simple_exec',&                                                                   ! executable
        &.true.,&                                                                          ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "model,filter,mask,compute"  )           ! GUI                                                      
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call multivol_assign%add_input(UI_SRCH, pgrp,                              gui_submenu="model",  gui_advanced=.false.)
        call multivol_assign%add_input(UI_SRCH, nsample,                           gui_submenu="search", gui_advanced=.false.)
        call multivol_assign%add_input(UI_SRCH, update_frac,                       gui_submenu="search", gui_advanced=.true.)
        call multivol_assign%add_input(UI_SRCH, nstates, required_override=.true., gui_submenu="search", gui_advanced=.false.)
        call multivol_assign%add_input(UI_SRCH, 'srch_oris', 'multi', 'Search orientations',&
        &'Search orientations(yes|no){yes}', '(yes|no){yes}', .false., 'single', gui_submenu="search", gui_advanced=.true.)
        ! filter controls
        call multivol_assign%add_input(UI_FILT, hp, gui_submenu="filter")
        call multivol_assign%add_input(UI_FILT, 'lpstart',  'num', 'Starting low-pass limit', 'Starting low-pass limit',&
            &'low-pass limit for the initial stage in Angstroms',  .false., 20., gui_submenu="filter")
        call multivol_assign%add_input(UI_FILT, 'lpstop',  'num', 'Final low-pass limit', 'Final low-pass limit',&
            &'low-pass limit for the final stage in Angstroms',    .false., 8., gui_submenu="filter")
        ! mask controls
        call multivol_assign%add_input(UI_MASK, mskdiam, gui_submenu="mask",    gui_advanced=.false.)
        ! computer controls
        call multivol_assign%add_input(UI_COMP, nparts,  required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call multivol_assign%add_input(UI_COMP, nthr,    gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('multivol_assign', multivol_assign, prgtab)
    end subroutine new_multivol_assign


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

end module simple_ui_api_abinitio3D
