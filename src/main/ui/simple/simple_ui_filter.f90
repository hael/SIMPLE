!@descr: module defining the user interfaces for filtering programs in the simple_exec suite
module simple_ui_filter
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('filter', 'Filtering', 80)
type(ui_program), target :: filter
type(ui_program), target :: uniform_filter2D
type(ui_program), target :: uniform_filter3D
type(ui_program), target :: nu_filt3D

contains

    subroutine construct_filter_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_filter(prgtab)
        call new_uniform_filter2D(prgtab)
        call new_uniform_filter3D(prgtab)
        call new_nu_filt3D(prgtab)
    end subroutine construct_filter_programs

subroutine new_filter( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call filter%new(&
        &'filter',&                                   ! name
        &'Apply a selected filter to an image stack or volume',& ! summary
        &'is a program for filtering stack/volume',&  ! help
        &'simple_exec',&                              ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call filter%add_input(UI_IMG, outstk, &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_IMG, outvol, &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call filter%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call filter%add_input(UI_IMG, 'stk',  'file', 'Stack to filter',  'Stack of images to filter', 'e.g. stk.mrcs',     .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_IMG, 'vol1', 'file', 'Volume to filter', 'Volume to filter',          'e.g. vol.mrc file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_requirement('input', 'Input data', 'Supply either a volume or an image stack.', &
            [character(len=4) :: 'stk ', 'vol1'], max_selected=1)
        ! search controls
        ! <empty>
        ! filter controls
        call filter%add_input(UI_FILT, lp, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, hp, &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'phrand', 'binary', 'Phase randomization', 'Fouirer phase randomization by white noise substitution(yes|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'bfac', 'num', 'B-factor of Gaussian low-/high-pass filter','B-factor of Gaussian low-/high-pass filter in Angstroms^2', 'B-factor in Angstroms^2{0}', .false., 0., &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'winsz', 'num', 'Half-window size', 'Half-window size(in pixels)', 'winsz in pixels', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'width', 'num', 'Cosine low-pass filter falloff',&
        &'Number of cosine edge pixels of Fourier low-pass filter in pixels', '# pixels cosine edge', .false., 10., &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'real_filter', 'multi', 'Real-space filter',&
        &'Real-space filter(median|average|stdev|bman|NLmean|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=7) :: 'median', 'average', 'stdev', 'bman', 'NLmean', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'fsc', 'file', 'FSC file', 'FSC file',          'e.g. fsc_state01.bin file', .false., '', &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, frcs, &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'filter', 'multi', 'Filter type(bs|nlmean|no){no}', 'Filter type(bs|nlmean|corr|no){no}','', .false., 'no', &
        &choices=ui_choices([character(len=6) :: 'bs', 'nlmean', 'no']), &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'lambda', 'num', 'BS smoother lambda', 'Strength of noise reduction', '(0.5-3.0){1.0}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        call filter%add_input(UI_FILT, 'sigma', 'num', 'sigma, for Gaussian generation', 'sigma, for Gaussian generation(in pixels)', &
        & '{1.}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        ! <empty>
        ! computer controls
        call filter%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('filter', filter, prgtab, UI_CATEGORY)
    end subroutine new_filter

    subroutine new_uniform_filter2D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_filter2D%new(&
        &'uniform_filter2D',&            ! name
        &'Filter 2D images using cross-validation',& ! summary
        &'is a program for 2D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! help
        &'simple_exec',&                 ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter2D%add_input(UI_IMG, 'stk',  'file', 'Odd stack',  'Odd stack',  'stack_even.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call uniform_filter2D%add_input(UI_IMG, 'stk2', 'file', 'Even stack', 'Even stack', 'stack_odd.mrc file',  .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call uniform_filter2D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter2D%add_input(UI_FILT, 'lpstart', 'num', 'Starting resolution limit', 'Starting resolution limit (in Angstroms)', 'in Angstroms', .true., -1., &
        &visibility=UI_VIS_STANDARD)
        call uniform_filter2D%add_input(UI_FILT, 'lpstop',  'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .true., -1., &
        &visibility=UI_VIS_STANDARD)
        ! mask controls
        call uniform_filter2D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call uniform_filter2D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('uniform_filter2D', uniform_filter2D, prgtab, UI_CATEGORY)
    end subroutine new_uniform_filter2D

    subroutine new_uniform_filter3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call uniform_filter3D%new(&
        &'uniform_filter3D',&                                   ! name
        &'Apply a uniform Butterworth filter to a 3D volume',& ! summary
        &'is a program for 3D uniform filter by minimizing/searching the fourier index of the CV cost function',& ! help
        &'simple_exec',&                                        ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call uniform_filter3D%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call uniform_filter3D%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call uniform_filter3D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call uniform_filter3D%add_input(UI_FILT, 'lpstart', 'num', 'Starting resolution limit', 'Starting resolution limit (in Angstroms)', 'in Angstroms', .true., -1., &
        &visibility=UI_VIS_STANDARD)
        call uniform_filter3D%add_input(UI_FILT, 'lpstop',  'num', 'Stopping resolution limit', 'Stopping resolution limit (in Angstroms)', 'in Angstroms', .true., -1., &
        &visibility=UI_VIS_STANDARD)
        call uniform_filter3D%add_input(UI_FILT, icm, &
        &visibility=UI_VIS_ADVANCED)
        call uniform_filter3D%add_input(UI_FILT, 'lambda', 'num', 'ICM lambda regularization parameter', 'Strength of noise reduction', '(0.01-3.0){1.0}', .false., 1.0, &
        &visibility=UI_VIS_ADVANCED)
        ! mask controls
        call uniform_filter3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        ! computer controls
        call uniform_filter3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('uniform_filter3D', uniform_filter3D, prgtab, UI_CATEGORY)
    end subroutine new_uniform_filter3D

    subroutine new_nu_filt3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call nu_filt3D%new(&
        &'nu_filt3D',&                                 ! name
        &'Apply spatially varying low-pass filtering to a 3D volume',& ! summary
        &'is a program for 3D nonuniform local low-pass filtering of even/odd volumes',& ! help
        &'simple_exec',&                                       ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call nu_filt3D%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call nu_filt3D%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call nu_filt3D%add_input(UI_IMG, outvol, required_override=.false., &
        &visibility=UI_VIS_ADVANCED)
        ! parameter input/output
        call nu_filt3D%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        call nu_filt3D%add_input(UI_MASK, mskdiam, &
        &visibility=UI_VIS_STANDARD)
        call nu_filt3D%add_input(UI_MASK, nu_envmsk, &
        &visibility=UI_VIS_ADVANCED)
        call nu_filt3D%add_input(UI_MASK, nu_msk_sig, &
        &visibility=UI_VIS_ADVANCED)
        ! nu_filt3D is the envelope-tuning program: these three are deliberately
        ! not offered anywhere else, where the NU_ENVMASK_* constants govern.
        call nu_filt3D%add_input(UI_MASK, nu_msk_beta, &
        &visibility=UI_VIS_ADVANCED)
        call nu_filt3D%add_input(UI_MASK, nu_msk_dens, &
        &visibility=UI_VIS_ADVANCED)
        call nu_filt3D%add_input(UI_MASK, nu_msk_rel, &
        &visibility=UI_VIS_ADVANCED)
        call nu_filt3D%add_input(UI_FILT, 'amsklp', 'num', 'NU envelope evidence scale',&
        &'Physical scale for smoothing the NU evidence margin, in Angstroms{8}', &
        &'evidence scale in Angstroms{8}', .false., 8., &
        &visibility=UI_VIS_ADVANCED)
        ! computer controls
        call nu_filt3D%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('nu_filt3D', nu_filt3D, prgtab, UI_CATEGORY)
    end subroutine new_nu_filt3D

end module simple_ui_filter
