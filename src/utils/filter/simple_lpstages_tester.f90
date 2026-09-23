!@descr: unit test routines for the low-pass and cropping schedules (mskdiam2lplimits, lpstages, lpstages_fast, lpstages_setlims) and the Butterworth kernel
! The clamps that turn a mask diameter into the 2D low-pass limits, the FRC-driven multi-stage schedule
! of refinement (stages from the FRC crossings, the linear fallback when the FRC never crosses, the crop
! box / sampling / shift-limit bookkeeping through the magic boxes), its two linear cousins, and the
! order-8 Butterworth transfer function against its closed form. References: lpstages_ref.py (stats
! batch scratch), a double-precision emulation of the same rules.
module simple_lpstages_tester
use simple_test_utils        ! assertions etc.
use simple_defs              ! sp
use simple_type_defs,        only: lp_crop_inf
use simple_string_utils,     only: int2str
use simple_estimate_ssnr,    only: mskdiam2lplimits, lpstages, lpstages_fast, lpstages_setlims
use simple_butterworth,      only: butterworth_filter
use simple_image,            only: image
implicit none
private
public :: run_all_lpstages_tests

integer, parameter :: BOX = 256, FILTSZ = BOX/2
real,    parameter :: SMPD = 1.3, LPSTART_LB = 10., LPSTART_DEFAULT = 20., LPFINAL = 6.
real,    parameter :: TOL = 1.e-4
! the falling FRC 1/(1+(k/40)^4): 5 stages from the crossings (lpstages_ref.py)
real,    parameter :: LP_K40(5)      = [20.0, 11.475862, 9.244444, 7.739535, 6.0]
real,    parameter :: CRIT_K40(5)    = [0.968405, 0.780955, 0.593505, 0.406055, 0.218605]
integer, parameter :: BOXCROP_K40(5) = [88, 108, 128, 150, 168]
! the faster-falling FRC 1/(1+(k/20)^4): the FRC value at the start limit (0.657) is above the particle
! threshold (0.65) but below the class-average one (0.80), so the two branches differ
real,    parameter :: LP_K20_PTCLS(4) = [20.0, 15.847619, 12.8, 6.0]
real,    parameter :: LP_K20_CAVGS(4) = [23.771429, 17.515789, 13.866667, 6.0]
real,    parameter :: CRIT_K20_CAVGS(4) = [0.80, 0.55, 0.30, 0.05]
integer, parameter :: BOXCROP_K20(4) = [88, 112, 140, 168]
! the linear fallback for a flat FRC, 10 stages from the start default to the final limit
integer, parameter :: BOXCROP_FLAT(10) = [88, 100, 108, 112, 128, 132, 140, 150, 160, 168]

contains

    subroutine run_all_lpstages_tests()
        write(*,'(A)') '**** running all low-pass stage tests ****'
        call test_mskdiam2lplimits()
        call test_lpstages_single_stage()
        call test_lpstages_from_frc()
        call test_lpstages_cavgs_vs_ptcls()
        call test_lpstages_flat_frc_fallback()
        call test_lpstages_fast()
        call test_lpstages_setlims()
        call test_butterworth()
    end subroutine run_all_lpstages_tests

    pure function frc_curve( k0 ) result( frc )
        real, intent(in) :: k0
        real    :: frc(FILTSZ)
        integer :: k
        do k = 1,FILTSZ
            frc(k) = 1. / (1. + (real(k)/k0)**4)
        end do
    end function frc_curve

    ! the crop bookkeeping every schedule shares: scale = box_crop/box, smpd_crop = smpd/scale, the shift
    ! limit is an alpha-helix width in cropped pixels clamped to [2,8], autoscale flags a real crop
    subroutine assert_crop_consistent( info, label )
        type(lp_crop_inf), intent(in) :: info
        character(len=*),  intent(in) :: label
        call assert_real(real(info%box_crop)/real(BOX), info%scale,     1.e-6, label//': scale = box_crop/box')
        call assert_real(SMPD/info%scale,               info%smpd_crop, 1.e-5, label//': smpd_crop = smpd/scale')
        call assert_real(min(8., max(2., 12.0/info%smpd_crop)), info%trslim, 1.e-5, label//': trslim from the helix width')
        call assert_true(info%l_autoscale .eqv. (info%box_crop < BOX), label//': l_autoscale flags a cropped box')
        call assert_true(info%l_lpset, label//': the low-pass limit is set')
    end subroutine assert_crop_consistent

    !---------------- mask diameter to 2D limits ----------------

    ! lpstart = clamp(d/12, 8, 15), lpstop = clamp(d/22, 5, 8), lpcen = clamp(d/6, 20, 30)
    subroutine test_mskdiam2lplimits()
        real :: lpstart, lpstop, lpcen
        write(*,'(A)') 'test_mskdiam2lplimits'
        call mskdiam2lplimits(300., lpstart, lpstop, lpcen)
        call assert_real(15., lpstart, TOL, '300 A: lpstart capped at 15')
        call assert_real(8.,  lpstop,  TOL, '300 A: lpstop capped at 8')
        call assert_real(30., lpcen,   TOL, '300 A: lpcen capped at 30')
        call mskdiam2lplimits(150., lpstart, lpstop, lpcen)
        call assert_real(12.5,     lpstart, TOL, '150 A: lpstart = d/12')
        call assert_real(6.818182, lpstop,  TOL, '150 A: lpstop = d/22')
        call assert_real(25.,      lpcen,   TOL, '150 A: lpcen = d/6')
        call mskdiam2lplimits(100., lpstart, lpstop, lpcen)
        call assert_real(8.333333, lpstart, TOL, '100 A: lpstart = d/12')
        call assert_real(5.,       lpstop,  TOL, '100 A: lpstop floored at 5')
        call assert_real(20.,      lpcen,   TOL, '100 A: lpcen floored at 20')
        call mskdiam2lplimits(60., lpstart, lpstop, lpcen)
        call assert_real(8.,  lpstart, TOL, '60 A: lpstart floored at 8')
        call assert_real(5.,  lpstop,  TOL, '60 A: lpstop floored at 5')
        call assert_real(20., lpcen,   TOL, '60 A: lpcen floored at 20')
    end subroutine test_mskdiam2lplimits

    !---------------- lpstages ----------------

    subroutine test_lpstages_single_stage()
        type(lp_crop_inf) :: info(1)
        real :: frc(FILTSZ)
        write(*,'(A)') 'test_lpstages_single_stage'
        frc = frc_curve(40.)
        call lpstages(BOX, 1, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, info, l_cavgs=.false., verbose=.false.)
        call assert_real(20., info(1)%lp, TOL, 'one stage: the start default (above the final limit)')
        call assert_int(88, info(1)%box_crop, 'one stage: 20 A targets smpd 6.7, the crop box is the 88 minimum')
        call assert_crop_consistent(info(1), 'one stage')
        call lpstages(BOX, 1, frc, SMPD, LPSTART_LB, 4., LPFINAL, info, l_cavgs=.false., verbose=.false.)
        call assert_real(20., info(1)%lp, TOL, 'one stage: a start default below the 20 A floor is floored')
        call lpstages(BOX, 1, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, 25., info, l_cavgs=.false., verbose=.false.)
        call assert_real(25., info(1)%lp, TOL, 'one stage: a final limit above the start wins')
    end subroutine test_lpstages_single_stage

    ! five stages from the FRC crossings of thresholds spaced evenly between the FRC at the start limit
    ! and at the final limit; the first stage is floored at 20 A, the last is the final limit; the crop
    ! boxes interpolate linearly through the magic boxes between the first and the last stage
    subroutine test_lpstages_from_frc()
        integer, parameter :: NSTAGES = 5
        type(lp_crop_inf) :: info(NSTAGES)
        real    :: frc(FILTSZ)
        integer :: i
        write(*,'(A)') 'test_lpstages_from_frc'
        frc = frc_curve(40.)
        call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, info, l_cavgs=.false., verbose=.false.)
        do i = 1,NSTAGES
            call assert_real(LP_K40(i),   info(i)%lp,       1.e-3, 'stage '//int2str(i)//': low-pass limit from the FRC crossing')
            call assert_real(CRIT_K40(i), info(i)%frc_crit, 1.e-4, 'stage '//int2str(i)//': critical FRC')
            call assert_int(BOXCROP_K40(i), info(i)%box_crop,     'stage '//int2str(i)//': crop box')
            call assert_crop_consistent(info(i), 'stage '//int2str(i))
        end do
        do i = 2,NSTAGES
            call assert_true(info(i)%lp < info(i-1)%lp, 'the limits tighten from stage '//int2str(i-1)//' to '//int2str(i))
            call assert_true(info(i)%box_crop >= info(i-1)%box_crop, 'the crop box grows from stage '//int2str(i-1)//' to '//int2str(i))
        end do
        call assert_real(3.781818, info(1)%smpd_crop, 1.e-4, 'stage 1: cropped sampling 1.3 * 256/88')
        call assert_real(3.173077, info(1)%trslim,    1.e-4, 'stage 1: shift limit 12 A in cropped pixels')
        call assert_real(6.057692, info(5)%trslim,    1.e-4, 'stage 5: shift limit at 1.98 A per pixel')
    end subroutine test_lpstages_from_frc

    ! the thresholds are floored at 0.65/0.03 for particles and 0.80/0.05 for class averages; with an FRC
    ! of 0.657 at the start limit the particle branch keeps the FRC value, the cavg branch takes 0.80
    subroutine test_lpstages_cavgs_vs_ptcls()
        integer, parameter :: NSTAGES = 4
        type(lp_crop_inf) :: info(NSTAGES)
        real    :: frc(FILTSZ)
        integer :: i
        write(*,'(A)') 'test_lpstages_cavgs_vs_ptcls'
        frc = frc_curve(20.)
        call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, info, l_cavgs=.false., verbose=.false.)
        do i = 1,NSTAGES
            call assert_real(LP_K20_PTCLS(i), info(i)%lp, 1.e-3, 'particles, stage '//int2str(i)//': low-pass limit')
            call assert_int(BOXCROP_K20(i), info(i)%box_crop, 'particles, stage '//int2str(i)//': crop box')
        end do
        call assert_real(0.657034, info(1)%frc_crit, 1.e-4, 'particles: the first threshold is the FRC at the start limit (above 0.65)')
        call assert_real(0.03,     info(4)%frc_crit, 1.e-6, 'particles: the last threshold is floored at 0.03')
        call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, info, l_cavgs=.true., verbose=.false.)
        do i = 1,NSTAGES
            call assert_real(LP_K20_CAVGS(i),   info(i)%lp,       1.e-3, 'cavgs, stage '//int2str(i)//': low-pass limit')
            call assert_real(CRIT_K20_CAVGS(i), info(i)%frc_crit, 1.e-6, 'cavgs, stage '//int2str(i)//': thresholds 0.80 to 0.05')
            call assert_int(BOXCROP_K20(i), info(i)%box_crop, 'cavgs, stage '//int2str(i)//': crop box')
        end do
        call assert_true(info(1)%lp > 20., 'cavgs: the stricter first threshold puts stage 1 above the floor')
    end subroutine test_lpstages_cavgs_vs_ptcls

    ! an FRC that never drops below a threshold gives no crossing: every stage is unset after the first pass
    ! and the schedule falls back to a linear ramp from the start default to the final limit
    subroutine test_lpstages_flat_frc_fallback()
        integer, parameter :: NSTAGES = 10
        type(lp_crop_inf) :: info(NSTAGES)
        real    :: frc(FILTSZ), expected
        integer :: i
        write(*,'(A)') 'test_lpstages_flat_frc_fallback'
        frc = 1.
        call lpstages(BOX, NSTAGES, frc, SMPD, LPSTART_LB, LPSTART_DEFAULT, LPFINAL, info, l_cavgs=.false., verbose=.false.)
        do i = 1,NSTAGES
            expected = LPSTART_DEFAULT - real(i-1) * (LPSTART_DEFAULT - LPFINAL) / real(NSTAGES-1)
            call assert_real(expected, info(i)%lp, 1.e-4, 'flat FRC, stage '//int2str(i)//': linear ramp')
            call assert_int(BOXCROP_FLAT(i), info(i)%box_crop, 'flat FRC, stage '//int2str(i)//': crop box')
            call assert_crop_consistent(info(i), 'flat FRC, stage '//int2str(i))
        end do
        call assert_real(1., info(3)%frc_crit, 1.e-6, 'flat FRC: the thresholds are the FRC value itself')
        ! a start default below the floor is floored before the ramp
        call lpstages(BOX, 3, frc, SMPD, LPSTART_LB, 12., LPFINAL, info(1:3), l_cavgs=.false., verbose=.false.)
        call assert_real(20., info(1)%lp, TOL, 'flat FRC: the ramp starts at the 20 A floor, not at 12')
        call assert_real(13., info(2)%lp, TOL, 'flat FRC: the ramp midpoint')
        call assert_real(6.,  info(3)%lp, TOL, 'flat FRC: the ramp ends at the final limit')
    end subroutine test_lpstages_flat_frc_fallback

    !---------------- lpstages_fast ----------------

    ! a linear ramp from max(lpstart, 20 A) to lpstop, each stage cropped on its own (smpd target 0.4 lp,
    ! at least 2.5 A), force_lpstart bypassing the floor
    subroutine test_lpstages_fast()
        type(lp_crop_inf) :: info(4)
        write(*,'(A)') 'test_lpstages_fast'
        call lpstages_fast(BOX, 4, SMPD, 30., 8., info, verbose=.false.)
        call assert_real(30.,       info(1)%lp, TOL, 'fast: starts at lpstart')
        call assert_real(22.666667, info(2)%lp, TOL, 'fast: linear ramp (2)')
        call assert_real(15.333333, info(3)%lp, TOL, 'fast: linear ramp (3)')
        call assert_real(8.,        info(4)%lp, TOL, 'fast: ends at lpstop')
        call assert_int(88,  info(1)%box_crop, 'fast: 30 A targets smpd 12, the crop box is the 88 minimum')
        call assert_int(88,  info(3)%box_crop, 'fast: 15.3 A targets smpd 6.1, still the minimum box')
        call assert_int(104, info(4)%box_crop, 'fast: 8 A targets smpd 3.2, magic box 104')
        call assert_real(3.2,  info(4)%smpd_crop, 1.e-5, 'fast: cropped sampling of the last stage')
        call assert_real(3.75, info(4)%trslim,    1.e-5, 'fast: shift limit of the last stage')
        call assert_crop_consistent(info(1), 'fast, stage 1')
        call assert_crop_consistent(info(4), 'fast, stage 4')
        call lpstages_fast(BOX, 3, SMPD, 12., 6., info(1:3), verbose=.false.)
        call assert_real(20., info(1)%lp, TOL, 'fast: lpstart below the floor is floored')
        call assert_real(13., info(2)%lp, TOL, 'fast: ramp from the floor')
        call lpstages_fast(BOX, 3, SMPD, 12., 6., info(1:3), force_lpstart=.true., verbose=.false.)
        call assert_real(12., info(1)%lp, TOL, 'fast: force_lpstart bypasses the floor')
        call assert_real(9.,  info(2)%lp, TOL, 'fast: ramp from the forced start')
        call lpstages_fast(BOX, 1, SMPD, 12., 6., info(1:1), verbose=.false.)
        call assert_real(20., info(1)%lp, TOL, 'fast, one stage: max(lpstop, floored lpstart)')
        call lpstages_fast(BOX, 1, SMPD, 12., 6., info(1:1), force_lpstart=.true., verbose=.false.)
        call assert_real(12., info(1)%lp, TOL, 'fast, one stage, forced: max(lpstop, lpstart)')
    end subroutine test_lpstages_fast

    !---------------- lpstages_setlims ----------------

    ! a linear ramp from lpstart to lpstop with no floor, crop boxes interpolated through the magic boxes
    subroutine test_lpstages_setlims()
        type(lp_crop_inf) :: info(5)
        integer, parameter :: BOXCROP(5) = [88, 108, 128, 150, 168]
        real,    parameter :: LPS(5)     = [15., 12.5, 10., 7.5, 5.]
        integer :: i
        write(*,'(A)') 'test_lpstages_setlims'
        call lpstages_setlims(BOX, 5, SMPD, 15., 5., info, verbose=.false.)
        do i = 1,5
            call assert_real(LPS(i), info(i)%lp, TOL, 'setlims, stage '//int2str(i)//': linear ramp, no floor')
            call assert_int(BOXCROP(i), info(i)%box_crop, 'setlims, stage '//int2str(i)//': crop box')
            call assert_crop_consistent(info(i), 'setlims, stage '//int2str(i))
        end do
        call lpstages_setlims(BOX, 1, SMPD, 15., 5., info(1:1), verbose=.false.)
        call assert_real(5.,   info(1)%lp, TOL, 'setlims, one stage: lpstop')
        call assert_int(168, info(1)%box_crop, 'setlims, one stage: cropped for lpstop')
        ! no cropping when the data are already coarser than the target sampling
        call lpstages_setlims(100, 1, 3.0, 15., 6., info(1:1), verbose=.false.)
        call assert_int(100, info(1)%box_crop, 'coarse data: the box is kept')
        call assert_real(1.,  info(1)%scale,     1.e-6, 'coarse data: scale 1')
        call assert_real(3.0, info(1)%smpd_crop, 1.e-6, 'coarse data: sampling kept')
        call assert_false(info(1)%l_autoscale, 'coarse data: no autoscaling')
    end subroutine test_lpstages_setlims

    !---------------- Butterworth ----------------

    ! |1/B_8(j s/fc)| = 1/sqrt(1 + (s/fc)^16): unity in the pass band, 1/sqrt(2) at the cut-off, the
    ! notch form is the product of a low-pass at c2 and a high-pass at c1
    subroutine test_butterworth()
        real, parameter :: EXPECTED(5) = [1.0, 0.999992, 0.707107, 0.165460, 0.003906] ! s = 0, 10, 20, 25, 40 for fc = 20
        real        :: fil(64), fil2(64), fil_rng(10:30)
        type(image) :: img
        integer     :: k
        logical     :: l_mono
        write(*,'(A)') 'test_butterworth'
        call butterworth_filter(20, fil)
        call assert_real(EXPECTED(1), fil(1),  1.e-6, 'Butterworth: unity at zero frequency')
        call assert_real(EXPECTED(2), fil(11), 1.e-4, 'Butterworth: flat pass band (s = fc/2)')
        call assert_real(EXPECTED(3), fil(21), 1.e-3, 'Butterworth: 1/sqrt(2) at the cut-off')
        call assert_real(EXPECTED(4), fil(26), 1.e-3, 'Butterworth: steep roll-off (s = 1.25 fc)')
        call assert_real(EXPECTED(5), fil(41), 1.e-4, 'Butterworth: stop band (s = 2 fc)')
        l_mono = .true.
        do k = 2,64
            if( fil(k) > fil(k-1) + 1.e-6 ) l_mono = .false.
        end do
        call assert_true(l_mono, 'Butterworth: monotone non-increasing')
        call assert_true(all(fil >= 0.) .and. all(fil <= 1.), 'Butterworth: values in [0,1]')
        call butterworth_filter(10, 30, fil2)
        call assert_real(0.995336, fil2(21), 1.e-3, 'notch 10-30: pass at s = 20')
        call assert_true(fil2(6) < 1.e-4, 'notch 10-30: blocked at s = 5')
        call assert_real(0.0996, fil2(41), 1.e-3, 'notch 10-30: 1/sqrt(1+(4/3)^16) at s = 40')
        call butterworth_filter(20, [10,30], fil_rng)
        call assert_real(fil(10), fil_rng(10), 1.e-6, 'ranged form: the same kernel over an index window (10)')
        call assert_real(fil(20), fil_rng(20), 1.e-6, 'ranged form: the same kernel over an index window (20)')
        call assert_real(fil(30), fil_rng(30), 1.e-6, 'ranged form: the same kernel over an index window (30)')
        ! the image form fills the same kernel and applies it
        call img%new([64,64,1], 1.)
        call img%set_rmat_at(33, 33, 1, 1.)
        call img%fft
        fil2 = 0.
        call butterworth_filter(img, 20, fil2)
        call assert_real(fil(21), fil2(21), 1.e-6, 'image form: the kernel it applied is the plain one')
        call assert_true(img%is_ft(), 'image form: the image stays in Fourier space')
        call img%kill
    end subroutine test_butterworth

end module simple_lpstages_tester
