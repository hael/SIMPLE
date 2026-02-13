!@descr: "filter" UI api (concrete implementation)
module simple_ui_api_filter
use simple_ui_api_modules
implicit none

type(ui_program), target :: filter
type(ui_program), target :: uniform_filter2D
type(ui_program), target :: uniform_filter3D

contains

    subroutine new_filter( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call filter%new(&
        &'filter',&                                   ! name
        &'Filter stack/volume',&                      ! descr_short
        &'is a program for filtering stack/volume',&  ! descr_long
        &'simple_exec',&                              ! executable
        &.false.)                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call filter%add_input(UI_IMG, outstk)
        call filter%add_input(UI_IMG, outvol)
        ! parameter input/output
        call filter%add_input(UI_PARM, smpd)
        ! alternative inputs
        call filter%add_input(UI_ALT, 'stk',  'file', 'Stack to filter',  'Stack of images to filter', 'e.g. stk.mrcs',     .false., '')
        call filter%add_input(UI_ALT, 'vol1', 'file', 'Volume to filter', 'Volume to filter',          'e.g. vol.mrc file', .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        call filter%add_input(UI_FILT, lp, required_override=.false.)
        call filter%add_input(UI_FILT, hp)
        call filter%add_input(UI_FILT, 'phrand', 'binary', 'Phase randomization', 'Fouirer phase randomization by white noise substitution(yes|no){no}', '(yes|no){no}', .false., 'no')
        call filter%add_input(UI_FILT, 'bfac', 'num', 'B-factor of Gaussian low-/high-pass filter','B-factor of Gaussian low-/high-pass filter in Angstroms^2', 'B-factor in Angstroms^2{0}', .false., 0.)
        call filter%add_input(UI_FILT, 'winsz', 'num', 'Half-window size', 'Half-window size(in pixels)', 'winsz in pixels', .false., 1.0)
        call filter%add_input(UI_FILT, 'width', 'num', 'Cosine low-pass filter falloff',&
        &'Number of cosine edge pixels of Fourier low-pass filter in pixels', '# pixels cosine edge', .false., 10.)
        call filter%add_input(UI_FILT, 'real_filter', 'multi', 'Real-space filter',&
        &'Real-space filter(median|average|stdev|bman|NLmean|no){no}', '(median|average|stdev|bman|NLmean|no){no}', .false., 'no')
        call filter%add_input(UI_FILT, 'fsc', 'file', 'FSC file', 'FSC file',          'e.g. fsc_state01.bin file', .false., '')
        call filter%add_input(UI_FILT, frcs)
        call filter%add_input(UI_FILT, 'filter', 'multi', 'Filter type(tv|nlmean|no){no}', 'Filter type(tv|nlmean|corr|no){no}', '(tv|nlmean|no){no}', .false., 'no')
        call filter%add_input(UI_FILT, 'lambda', 'num', 'Tv filter lambda', 'Strength of noise reduction', '(0.5-3.0){1.0}', .false., 1.0)
        call filter%add_input(UI_FILT, 'sigma', 'num', 'sigma, for Gaussian generation', 'sigma, for Gaussian generation(in pixels)', &
        & '{1.}', .false., 1.0)
        ! mask controls
        ! <empty>
        ! computer controls
        call filter%add_input(UI_COMP, nthr)
        call add_ui_program('filter', filter, prgtab)
    end subroutine new_filter


    subroutine new_uniform_filter2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_filter2D%new(&
        &'uniform_filter2D',&            ! name
        &'Uniform 2D filter',&           ! descr_short
        &'is a program for 2D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                 ! executable
        &.false.)                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter2D%add_input(UI_IMG, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '')
        call uniform_filter2D%add_input(UI_IMG, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '')
        ! parameter input/output
        call uniform_filter2D%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter2D%add_input(UI_FILT, 'lpstart', 'num', 'Starting resolution limit', 'Starting resolution limit (in Angstroms)', 'in Angstroms', .true., -1.)
        call uniform_filter2D%add_input(UI_FILT, 'lpstop',  'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .true., -1.)
        ! mask controls
        call uniform_filter2D%add_input(UI_MASK, mskdiam)
        ! computer controls
        call uniform_filter2D%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('uniform_filter2D', uniform_filter2D, prgtab)
    end subroutine new_uniform_filter2D

    subroutine new_uniform_filter3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_filter3D%new(&
        &'uniform_filter3D',&                                   ! name
        &'Uniform Butterworth 3D filter',&                      ! descr_short
        &'is a program for 3D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! descr_long
        &'simple_exec',&                                        ! executable
        &.false.)                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter3D%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call uniform_filter3D%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call uniform_filter3D%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter3D%add_input(UI_FILT, 'lpstart', 'num', 'Starting resolution limit', 'Starting resolution limit (in Angstroms)', 'in Angstroms', .true., -1.)
        call uniform_filter3D%add_input(UI_FILT, 'lpstop',  'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .true., -1.)
        call uniform_filter3D%add_input(UI_FILT, icm)
        call uniform_filter3D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0)
        ! mask controls
        call uniform_filter3D%add_input(UI_MASK, mskdiam)
        call uniform_filter3D%add_input(UI_MASK, mskfile)
        ! computer controls
        call uniform_filter3D%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('uniform_filter3D', uniform_filter3D, prgtab)
    end subroutine new_uniform_filter3D

end module simple_ui_api_filter
