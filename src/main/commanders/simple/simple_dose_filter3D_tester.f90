!@descr: integration of dose_filter3D with a reconstructed fraction dose-series manifest
module simple_dose_filter3D_tester
use simple_defs,                        only: CWD_GLOB
use simple_string,                      only: string
use simple_cmdline,                     only: cmdline
use simple_fileio,                      only: file_exists
use simple_syslib,                      only: simple_chdir
use simple_commanders_resolest,         only: commander_dose_filter3D
use simple_dose_curve,                  only: dose_curve
use simple_refine_motion_model_tester,  only: check_fraction_dose_export
use simple_flex_pca_application_tester, only: create_flex_pca_phantom_fixture, FLEX_PHANTOM_SMPD, FLEX_PHANTOM_MSKDIAM
use simple_test_utils,                  only: assert_true, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_dose_filter3D_tests

integer, parameter :: RECON_NPTCLS = 120

contains

    !> Run the dose consumer on a fresh producer fixture, independently of the producer test suite.
    subroutine run_all_dose_filter3D_tests()
        class(cmdline), allocatable :: filter_cline
        type(commander_dose_filter3D) :: filter
        type(dose_curve) :: curve
        type(string) :: cwd, root, arm, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, message
        integer :: nfail0
        nfail0 = tests_failed
        call enter_fixture('dose_filter3d', cwd, root)
        allocate(filter_cline)
        call create_flex_pca_phantom_fixture(root, 1, truth_a, truth_b, truth_mean, truth_diff, oritab, stack, nptcls=RECON_NPTCLS)
        call check_fraction_dose_export(root, 'valid', stack, oritab, .true., .false.)
        arm = root//'/valid'
        if( file_exists(arm//'/dose_series.txt') )then
            call simple_chdir(arm)
            CWD_GLOB = arm%to_char()
            call filter_cline%set('prg', 'dose_filter3D')
            call filter_cline%set('filetab', arm//'/dose_series.txt')
            call filter_cline%set('smpd', FLEX_PHANTOM_SMPD)
            call filter_cline%set('mskdiam', FLEX_PHANTOM_MSKDIAM)
            call filter_cline%set('hp', 30.)
            call filter_cline%set('lp', 8.)
            call filter_cline%set('nthr', 1)
            call filter_cline%set('mkdir', 'no')
            call filter_cline%set('outvol', arm//'/dose_result.mrc')
            call filter%execute(filter_cline)
            call assert_true(file_exists(arm//'/dose_result_even.mrc'), 'dose filter3D: produced even accumulation')
            call assert_true(file_exists(arm//'/dose_result_odd.mrc'), 'dose filter3D: produced odd accumulation')
            call assert_true(curve%read(arm//'/dose_result_curve.txt', message), 'dose filter3D: wrote a dose curve reextract reads')
            call assert_true(file_exists(arm//'/dose_result_fsc_masked.bin'), 'dose filter3D: wrote the masked FSC')
            call assert_true(file_exists(arm//'/dose_result_fsc.txt'), 'dose filter3D: wrote the FSC table')
            call curve%kill
        endif
        call filter_cline%kill
        call simple_chdir(root)
        CWD_GLOB = root%to_char()
        call leave_fixture(cwd, root, nfail0)
        CWD_GLOB = cwd%to_char()
    end subroutine run_all_dose_filter3D_tests

end module simple_dose_filter3D_tester
