!@descr: "single_ui_api_nano2D" UI api (concrete implementation)
module single_ui_api_nano2D
use simple_ui_api_modules
implicit none

type(ui_program), target :: analysis2D_nano
type(ui_program), target :: center2D_nano
type(ui_program), target :: cluster2D_nano
type(ui_program), target :: estimate_diam

contains

    subroutine register_single_ui_nano2D(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('analysis2D_nano', analysis2D_nano, prgtab)
        call add_ui_program('center2D_nano', center2D_nano, prgtab)
        call add_ui_program('cluster2D_nano', cluster2D_nano, prgtab)
        call add_ui_program('estimate_diam', estimate_diam, prgtab)
    end subroutine register_single_ui_nano2D

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_analysis2D_nano
        ! PROGRAM SPECIFICATION
        call analysis2D_nano%new(&
        &'analysis2D_nano', &                                         ! name
        &'2D analysis (centering, diameter estimation & clustering) for nanocrystal time-series',& ! descr_short
        &'is a program for 2D analysis for nanycrystal time-series',& ! descr long
        &'single_exec',&                                              ! executable
        &.true., gui_advanced=.false.)                                ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call analysis2D_nano%add_input(UI_PARM, element)
        ! alternative inputs
        ! <empty>
        ! search controls
        call analysis2D_nano%add_input(UI_SRCH, nptcls_per_cls)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call analysis2D_nano%add_input(UI_COMP, nthr)
        call analysis2D_nano%add_input(UI_COMP, script)
    end subroutine new_analysis2D_nano

    subroutine new_center2D_nano
        ! PROGRAM SPECIFICATION
        call center2D_nano%new(&
        &'center2D_nano',&                                                      ! name
        &'Simultaneous 2D alignment and clustering of nanoparticle images',&    ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm&
        & suitable for the first pass of cleanup after time-series tracking',&  ! descr_long
        &'single_exec',&                                                        ! executable
        &.true., gui_advanced=.false.)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call center2D_nano%add_input(UI_SRCH, ncls, required_override=.true.)
        call center2D_nano%add_input(UI_SRCH, trs)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call center2D_nano%add_input(UI_COMP, nthr)
        call center2D_nano%add_input(UI_COMP, script)
    end subroutine new_center2D_nano

    subroutine new_cluster2D_nano
        ! PROGRAM SPECIFICATION
        call cluster2D_nano%new(&
        &'cluster2D_nano',&                                                                 ! name
        &'Simultaneous 2D alignment and clustering of time-series of nanoparticle images',& ! descr_short
        &'is a distributed workflow implementing a reference-free 2D alignment/clustering algorithm for time-series of nanoparticle images',& ! descr_long
        &'single_exec',&                                                                    ! executable
        &.true., gui_advanced=.false.)                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call cluster2D_nano%add_input(UI_PARM, moldiam)
        ! alternative inputs
        ! <empty>
        ! search controls
        call cluster2D_nano%add_input(UI_SRCH, nptcls_per_cls)
        call cluster2D_nano%add_input(UI_SRCH, 'center', 'binary', 'Center class averages', 'Center class averages by their center of &
        &gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call cluster2D_nano%add_input(UI_SRCH, 'winsz', 'num', 'Half-window size', 'Half-window size(frames)', 'winsz in # frames', .false., 3.0)
        call cluster2D_nano%add_input(UI_SRCH, maxits)
        call cluster2D_nano%add_input(UI_SRCH, trs)
        ! filter controls
        call cluster2D_nano%add_input(UI_FILT, hp)
        call cluster2D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the class averages and centering', 'centering low-pass limit in &
        &Angstroms{5.0}', .false., 5.)
        call cluster2D_nano%add_input(UI_FILT, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit{1.0}', 'low-pass limit in Angstroms', .false., 1.)
        call cluster2D_nano%add_input(UI_FILT, 'lpstart', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'initial low-pass limit in Angstroms', .false., 1.)
        call cluster2D_nano%add_input(UI_FILT, 'lpstop', 'num', 'Final low-pass limit', 'Final low-pass limit{1.0}', 'final low-pass limit in Angstroms', .false., 1.)
        ! mask controls
        call cluster2D_nano%add_input(UI_MASK, mskdiam)
        ! computer controls
        call cluster2D_nano%add_input(UI_COMP, nparts, required_override=.false.)
        call cluster2D_nano%add_input(UI_COMP, nthr)
        call cluster2D_nano%add_input(UI_COMP, script)
    end subroutine new_cluster2D_nano

    subroutine new_estimate_diam
        ! PROGRAM SPECIFICATION
        call estimate_diam%new(&
        &'estimate_diam',&                                                                                              ! name
        &'Estimation of a suitable mask diameter for nanoparticle time-series',&                                        ! descr_short
        &'is a program for estimation of a suitable mask diameter for spherical masking of nanoparticle time-series ',& ! descr_long
        &'single_exec',&                                                                                                ! executable
        &.false., gui_advanced=.false.)                                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call estimate_diam%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call estimate_diam%add_input(UI_PARM, smpd)
        call estimate_diam%add_input(UI_PARM, 'roavg', 'binary', 'Rotationally average', 'Rotationally average before diam estimate(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call estimate_diam%add_input(UI_FILT, lp, descr_short_override='low-pass limit in Angstroms{7.}')
        ! mask controls
        call estimate_diam%add_input(UI_MASK, mskdiam)
        ! computer controls
        call estimate_diam%add_input(UI_COMP, nthr)
    end subroutine new_estimate_diam

end module single_ui_api_nano2D
