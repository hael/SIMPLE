!@descr: unit tests of the abinitio3D_addon validation report (simple_abinitio3D_addon_report)
! The verdict on synthetic logistic FSC curves (improved and regressed beyond
! one FSC=0.143 shell, unchanged within it, curves of another length not
! compared); a map against itself and its cohort stand-in (correlation one,
! no docking, cohort resolutions found); and the file round trip of stage,
! state and cohort records. One small synthetic map; the docking branch takes
! seconds and is the library sub-suite of simple_abinitio3D_addon_dock_tester.
module simple_abinitio3D_addon_report_tester
use simple_core_module_api
use simple_image,                   only: image
use simple_abinitio3D_addon_report, only: abinitio3D_addon_report
use simple_test_utils
implicit none

public :: run_all_abinitio3D_addon_report_tests
private
#include "simple_local_flags.inc"

real,             parameter :: SMPD     = 2.
real,             parameter :: MSKDIAM  = 100.
integer,          parameter :: BOX_FSC  = 64
integer,          parameter :: BOX_MAP  = 64
integer,          parameter :: SEED     = 20260927
character(len=*), parameter :: F_MAP    = 'tmp_addon_report_map.mrc'
character(len=*), parameter :: F_REPORT = 'tmp_addon_report.txt'

contains

    subroutine run_all_abinitio3D_addon_report_tests()
        write(logfhandle,'(A)') '**** running all abinitio3D add-on report tests ****'
        call write_fixtures()
        call test_fsc_verdicts()
        call test_map_comparison()
        call test_round_trip()
        call del_file(F_MAP)
        call del_file(F_REPORT)
    end subroutine run_all_abinitio3D_addon_report_tests

    ! ---- fixtures -------------------------------------------------------------

    !> an asymmetric map of four blobs over weak noise
    subroutine write_fixtures()
        real, parameter :: CENTRES(3,4) = reshape([-10.,4.,2., 9.,-6.,0., 2.,11.,-8., -3.,-9.,7.], [3,4])
        real, allocatable :: rmat(:,:,:), noise(:,:,:)
        real    :: d2
        integer :: i, j, k, b, c
        call set_fixed_seed(SEED)
        allocate(rmat(BOX_MAP,BOX_MAP,BOX_MAP), source=0.)
        allocate(noise(BOX_MAP,BOX_MAP,BOX_MAP))
        call random_number(noise)
        c = BOX_MAP/2 + 1
        do k = 1, BOX_MAP
            do j = 1, BOX_MAP
                do i = 1, BOX_MAP
                    do b = 1, 4
                        d2 = (real(i-c) - CENTRES(1,b))**2 + (real(j-c) - CENTRES(2,b))**2 + (real(k-c) - CENTRES(3,b))**2
                        rmat(i,j,k) = rmat(i,j,k) + real(b) * exp(-0.5 * d2 / 9.)
                    enddo
                enddo
            enddo
        enddo
        rmat = rmat + 0.02 * (noise - 0.5)
        call write_map(F_MAP, rmat)
    end subroutine write_fixtures

    subroutine write_map( fname, rmat )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: rmat(:,:,:)
        type(image) :: vol
        call vol%new([BOX_MAP,BOX_MAP,BOX_MAP], SMPD, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(fname), del_if_exists=.true.)
        call vol%kill
    end subroutine write_map

    !> a logistic FSC falling through 0.5 near shell k0
    function fsc_curve( n, k0 ) result( fsc )
        integer, intent(in) :: n
        real,    intent(in) :: k0
        real :: fsc(n)
        integer :: k
        fsc = [(1. / (1. + exp((real(k) - k0) / 1.5)), k = 1, n)]
    end function fsc_curve

    ! ---- tests ----------------------------------------------------------------

    subroutine test_fsc_verdicts()
        type(abinitio3D_addon_report) :: report
        real, allocatable :: res(:), base(:)
        integer :: n
        write(logfhandle,'(A)') '-- test_fsc_verdicts'
        res  = get_resarr(BOX_FSC, SMPD)
        n    = size(res)
        base = fsc_curve(n, 15.)
        call report%new(5, 8, SMPD, MSKDIAM)
        call report%compare_fsc(1, base, fsc_curve(n, 20.), res)
        call report%compare_fsc(2, base, fsc_curve(n, 15.), res)
        call report%compare_fsc(3, base, fsc_curve(n, 16.), res)
        call assert_char('IMPROVED',  trim(report%get_verdict(1)), 'five shells finer is an improvement')
        call assert_int(5, report%get_dshell(1), 'the FSC=0.143 shell moved by five')
        call assert_char('UNCHANGED', trim(report%get_verdict(2)), 'the same curve is unchanged')
        call assert_char('UNCHANGED', trim(report%get_verdict(3)), 'one shell is within the tolerance')
        call assert_false(report%any_regressed(), 'no state regressed yet')
        call report%compare_fsc(4, base, fsc_curve(n, 10.), res)
        call assert_char('REGRESSED', trim(report%get_verdict(4)), 'five shells coarser is a regression')
        call assert_true(report%any_regressed(), 'a regressed state is reported')
        call report%compare_fsc(5, base, fsc_curve(n - 1, 15.), res)
        call assert_char('NOT_COMPARED', trim(report%get_verdict(5)), 'curves of another length are not compared')
        call report%kill
    end subroutine test_fsc_verdicts

    subroutine test_map_comparison()
        type(abinitio3D_addon_report) :: report
        real, allocatable :: res(:)
        write(logfhandle,'(A)') '-- test_map_comparison'
        res = get_resarr(BOX_FSC, SMPD)
        call report%new(1, 8, SMPD, MSKDIAM)
        call report%compare_fsc(1, fsc_curve(size(res), 15.), fsc_curve(size(res), 18.), res)
        call report%compare_maps(1, string(F_MAP), string(F_MAP))
        call assert_real(1., report%get_corr(1), 1.e-4, 'a map correlates fully with itself')
        call assert_real(-1., report%get_dock_angle(1), 0., 'a map in the base frame is not docked')
        call report%compare_cohort(1, string(F_MAP), string(F_MAP))
        call assert_true(report%get_cohort_res0143(1) > 0., 'the cohort map resolution is found')
        call report%compare_maps(1, string(F_MAP), string('tmp_addon_report_absent.mrc'))
        call assert_real(0., report%get_corr(1), 0., 'a missing map is not compared')
        call report%kill
    end subroutine test_map_comparison

    subroutine test_round_trip()
        type(abinitio3D_addon_report) :: report, back
        real, allocatable :: res(:)
        write(logfhandle,'(A)') '-- test_round_trip'
        res = get_resarr(BOX_FSC, SMPD)
        call report%new(2, 8, SMPD, MSKDIAM)
        call report%set_stage(3, 12., 12., 12., 12.)
        call report%set_stage(4, 8.5, 8.5, 10., 10.)
        call report%set_stage(6, 6., 0., -1., -1.)
        call report%set_populations(1, 1200, 600)
        call report%set_populations(2, 900, 450)
        call report%compare_fsc(1, fsc_curve(size(res), 15.), fsc_curve(size(res), 20.), res)
        call report%compare_fsc(2, fsc_curve(size(res), 15.), fsc_curve(size(res), 10.), res)
        call report%compare_maps(1, string(F_MAP), string(F_MAP))
        call report%compare_cohort(1, string(F_MAP), string(F_MAP))
        call report%write(string(F_REPORT))
        call back%read(string(F_REPORT))
        call assert_char(trim(report%get_verdict(1)), trim(back%get_verdict(1)), 'state 1 verdict read back')
        call assert_char(trim(report%get_verdict(2)), trim(back%get_verdict(2)), 'state 2 verdict read back')
        call assert_int(report%get_dshell(2), back%get_dshell(2), 'shell move read back')
        call assert_real(report%get_corr(1), back%get_corr(1), 1.e-4, 'map correlation read back')
        call assert_real(report%get_cohort_res0143(1), back%get_cohort_res0143(1), 1.e-3, 'cohort resolution read back')
        call assert_true(back%any_regressed(), 'the regressed state is read back')
        call report%kill
        call back%kill
    end subroutine test_round_trip

end module simple_abinitio3D_addon_report_tester
