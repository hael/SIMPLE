!@descr: unit tests for the dose-series attenuation model, its dose curve, its accumulations and its cross-validation
! Synthetic series on a small shell grid: the model of doc/implementation_notes/planned/dose_filter3D.md
! is generated exactly, so every expectation here is a closed form, not a recorded run.
module simple_dose_series_tester
use simple_test_utils      ! assertions etc.
use simple_defs            ! sp, dp
use simple_string,      only: string
use simple_dose_series, only: dose_series, shell_independence, shell_correlation, line_tau2, DOSE_WEIGHTS_MODEL, DOSE_WEIGHTS_CONTROL, &
                              &DOSE_WEIGHTS_FLAT
use simple_dose_curve,  only: dose_curve
!$ use omp_lib,         only: omp_get_max_threads, omp_get_num_procs, omp_set_num_threads
implicit none
private
public :: run_all_dose_series_tests

integer,  parameter :: NSH = 24, NPERSH = 256
real(dp), parameter :: SMPD = 1.5_dp
integer,  parameter :: BOX  = 2*NSH
real(dp), parameter :: TOL_LOOSE = 0.15_dp

contains

    subroutine run_all_dose_series_tests()
        write(*,'(A)') '**** running all dose series tests ****'
        call test_noise_free_attenuation()
        call test_recovers_the_dose_curve()
        call test_no_damage_gives_precision_weights()
        call test_two_windows_are_identifiable()
        call test_accumulation_is_a_weighted_mean()
        call test_band_follows_the_best_window()
        call test_statistics_are_scale_free()
        call test_control_and_flat_weights()
        call test_zerodose_rescaling()
        call test_paired_cv_rejects_a_wrong_control()
        call test_paired_cv_keeps_a_true_control()
        call test_parallel_matches_serial()
        call test_errors_are_calibrated()
        call test_exposure_errors_share_the_reference()
        call test_relative_b_errors_share_the_reference()
        call test_shell_independence()
        call test_failed_fits_are_rejected()
        call test_line_tau2()
        call test_shell_exposure()
        call test_exports_the_dose_curve()
        call test_step_exposure()
        call test_reports_against_the_sharpest()
    end subroutine run_all_dose_series_tests

    !> A window sharper than window 1 becomes the reference of the model table: it reads c = 1 and B = 0, window 1
    !! reads the positive difference, no B is negative, and the dose slope is the fit's
    subroutine test_reports_against_the_sharpest()
        integer,          parameter :: NW = 5
        character(len=*), parameter :: FNAME = 'dose_series_tester_sharpest.txt'
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(NW), bw(NW), sigma(NW), rows(10,NW), line(4), beta
        character(len=512) :: text
        integer(kind=8) :: seed
        integer         :: w, funit, ios
        write(*,'(A)') 'test_reports_against_the_sharpest'
        call shell_grid(shell, kfreq)
        bw    = [0._dp, -6._dp, 10._dp, 30._dp, 50._dp]
        cw    = [1._dp, 1.02_dp, 0.98_dp, 0.95_dp, 0.92_dp]
        sigma = 0.01_dp
        seed  = 31415_8
        call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(NW, kfreq, 300._dp, .true.)
        do w = 1, NW
            call ds%set_window(w, real(2*w-2,dp), real(2*w,dp))
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        beta = ds%get_beta()
        call ds%write_model_table(string(FNAME))
        call ds%kill
        open(newunit=funit, file=FNAME, status='old', action='read', iostat=ios)
        call assert_true(ios == 0, 'the model table is written')
        if( ios /= 0 ) return
        read(funit,'(A)') text
        call assert_true(index(text, 'relative to window 2') > 0, 'the table names the sharpest window as its reference')
        read(funit,*)
        do w = 1, NW
            read(funit,*) rows(:,w)
        enddo
        read(funit,'(A)') text
        read(text(index(text,'tau')+3:),*) line
        close(funit, status='delete')
        ! the shift is the reference's B*, which its own B_w differs from by its shrinkage
        call assert_true(abs(rows(7,2)) < 1.e-9_dp .and. abs(rows(3,2) - 1._dp) < 1.e-9_dp .and. rows(5,2) > 0._dp, &
            &'the reference reads c = 1 and B* = 0, with uncertainty on its measured B residual')
        call assert_true(abs(rows(4,2)) < 0.5_dp, 'the reference window B_w is its shrinkage from zero')
        call assert_true(abs(rows(4,1) - 6._dp) < 1._dp, 'window 1 reads the B it lacks against the sharpest')
        call assert_true(all(rows(7,:) >= 0._dp), 'no B* is negative against the sharpest window')
        call assert_true(abs(line(2) - beta) < 1.e-6_dp, 'the dose slope is the fit, whatever the reference')
        deallocate(em, om, shell, kfreq)
    end subroutine test_reports_against_the_sharpest

    !> With B rising as the square of dose the fading is not exponential in dose: the critical exposure between
    !! neighbouring windows falls along the series as 2/(b2 k^2 (D_w + D_{w+1} - 2 D_1)). The per-window model
    !! read by the measured estimator still matches the measured critical exposure; the lines' slope does not.
    subroutine test_step_exposure()
        integer,          parameter :: NW = 5
        character(len=*), parameter :: STEPS_FNAME = 'dose_series_tester_steps.txt'
        character(len=*), parameter :: DECAY_FNAME = 'dose_series_tester_steps_decay.txt'
        real(dp),         parameter :: B2 = 1._dp, RES(2) = [12._dp, 6._dp]
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(NW), bw(NW), sigma(NW), dose(NW), row(2+3*size(RES)), vals(7), truth, err, err_est, gap
        integer(kind=8) :: seed
        integer         :: w, i, s, funit, ios, nfinite, nrow
        write(*,'(A)') 'test_step_exposure'
        call shell_grid(shell, kfreq)
        do w = 1, NW
            dose(w)  = real(2*w - 1, dp)
            bw(w)    = B2 * (dose(w) - dose(1))**2
            cw(w)    = 1._dp
            sigma(w) = 0.002_dp
        enddo
        seed = 5150_8
        call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(NW, kfreq, 300._dp, .true.)
        do w = 1, NW
            call ds%set_window(w, dose(w) - 1._dp, dose(w) + 1._dp)
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%write_step_exposure(string(STEPS_FNAME), RES)
        call ds%write_decay_table(string(DECAY_FNAME))
        call ds%kill
        open(newunit=funit, file=STEPS_FNAME, status='old', action='read', iostat=ios)
        call assert_true(ios == 0, 'the step exposure table is written')
        if( ios == 0 )then
            read(funit,*)
            err = 0._dp; nfinite = 0
            do w = 1, NW - 1
                read(funit,*,iostat=ios) row
                if( ios /= 0 ) exit
                do i = 1, size(RES)
                    if( row(3*i) /= row(3*i) ) cycle                ! NaN: not measured
                    s       = nint(1._dp / (RES(i) * kfreq(1)))
                    truth   = 2._dp / (B2 * kfreq(s)**2 * (dose(w) + dose(w+1) - 2._dp*dose(1)))
                    nfinite = nfinite + 1
                    err     = max(err, abs(row(3*i) - truth) / truth)
                enddo
            enddo
            close(funit, status='delete')
            write(logfhandle,'(A,ES12.3)') 'neighbouring-window critical exposure, worst relative error: ', err
            call assert_true(nfinite == size(RES)*(NW-1), 'every neighbouring pair is measured at both resolutions')
            call assert_true(err < 0.05_dp, 'the neighbouring-window critical exposure follows the dose dependence')
        endif
        open(newunit=funit, file=DECAY_FNAME, status='old', action='read', iostat=ios)
        call assert_true(ios == 0, 'the decay table is written')
        if( ios == 0 )then
            read(funit,*); read(funit,*)
            nrow = 0; err_est = 0._dp; gap = 0._dp
            do
                read(funit,*,iostat=ios) vals
                if( ios /= 0 ) exit
                if( vals(3) /= vals(3) .or. vals(7) /= vals(7) ) cycle
                nrow    = nrow + 1
                err_est = max(err_est, abs(vals(7) - vals(3)) / vals(3))
                gap     = max(gap,     abs(vals(5) - vals(3)) / vals(3))
            enddo
            close(funit, status='delete')
            write(logfhandle,'(A,2ES12.3)') 'accelerating fading, per-window model and lines against measured: ', err_est, gap
            call assert_true(nrow >= 6, 'the critical exposure is measured in most shells')
            call assert_true(err_est < 0.05_dp, 'the per-window model read by the same estimator matches the measured')
            call assert_true(gap > 0.2_dp, 'the lines measure another thing when the fading accelerates')
        endif
        deallocate(em, om, shell, kfreq)
    end subroutine test_step_exposure

    !> The exported curve holds each window's contrast and B-factor, and past the last window the lines of both
    !! in dose: read back per movie fraction from the file reextract reads. The decay table's model critical
    !! exposure fades by both lines, 1/(2 gamma + beta k^2/2).
    subroutine test_exports_the_dose_curve()
        integer,          parameter :: NW = 5, NFRAC = 14
        character(len=*), parameter :: FNAME = 'dose_series_tester_curve.txt'
        character(len=*), parameter :: DECAY_FNAME = 'dose_series_tester_curve_decay.txt'
        real(dp),         parameter :: BETA_TRUE = 4._dp, GAMMA_TRUE = 0.02_dp
        type(dose_series)        :: ds
        type(dose_curve)         :: curve, readback
        type(string)             :: message
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), c(:), b(:), c2(:), b2(:)
        real(dp)        :: cw(NW), bw(NW), sigma(NW), dose(NW), x, err_c, err_b, vals(7), truth, err_ne, err_est
        integer(kind=8) :: seed
        integer         :: w, f, ffrom(NW), fto(NW), funit, ios, nrow
        logical         :: ok
        write(*,'(A)') 'test_exports_the_dose_curve'
        call shell_grid(shell, kfreq)
        ! windows of two fractions of 1 e/A^2 at 300 kV, centred at 1, 3, 5, 7 and 9
        do w = 1, NW
            ffrom(w) = 2*w - 1
            fto(w)   = 2*w
            dose(w)  = real(2*w - 1, dp)
            bw(w)    = BETA_TRUE * (dose(w) - dose(1))
            cw(w)    = exp(-GAMMA_TRUE * (dose(w) - dose(1)))
            sigma(w) = 0.01_dp
        enddo
        seed = 4242_8
        call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(NW, kfreq, 300._dp, .true.)
        do w = 1, NW
            call ds%set_window(w, real(ffrom(w) - 1, dp), real(fto(w), dp))
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%export_curve(curve, ffrom, fto, 1._dp)
        call ds%write_decay_table(string(DECAY_FNAME))
        call ds%kill
        open(newunit=funit, file=DECAY_FNAME, status='old', action='read', iostat=ios)
        call assert_true(ios == 0, 'the decay table is written')
        if( ios == 0 )then
            read(funit,*); read(funit,*)
            nrow = 0; err_ne = 0._dp; err_est = 0._dp
            do
                read(funit,*,iostat=ios) vals
                if( ios /= 0 ) exit
                nrow    = nrow + 1
                truth   = 1._dp / (2._dp*GAMMA_TRUE + BETA_TRUE*vals(2)**2/2._dp)
                err_ne  = max(err_ne,  abs(vals(5) - truth) / truth)
                err_est = max(err_est, abs(vals(7) - truth) / truth)
            enddo
            close(funit, status='delete')
            write(logfhandle,'(A,2ES12.3)') 'model critical exposure, lines and per window, worst relative error: ', &
                &err_ne, err_est
            call assert_true(nrow >= 6, 'the decay table has a row per band shell')
            call assert_true(err_ne < 0.02_dp, 'the model critical exposure fades by the contrast and B-factor lines')
            call assert_true(err_est < 0.02_dp, 'on fading exponential in dose the per-window model gives the same')
        endif
        call curve%write(string(FNAME))
        ok = readback%read(string(FNAME), message)
        call assert_true(ok, 'the exported dose curve reads back')
        if( ok )then
            call curve%fraction_weights(NFRAC, 1._dp, 300._dp, c, b)
            call readback%fraction_weights(NFRAC, 1._dp, 300._dp, c2, b2)
            call assert_true(all(c == c2) .and. all(b == b2), 'the file holds the exported curve exactly')
            ! the truth is linear in dose, so the fractions between and past the windows recover it alike
            err_c = 0._dp; err_b = 0._dp
            do f = 1, NFRAC
                x     = real(f, dp) - 0.5_dp - dose(1)
                err_c = max(err_c, abs(log(c(f)) + GAMMA_TRUE * x))
                err_b = max(err_b, abs(b(f) - BETA_TRUE * x))
            enddo
            write(logfhandle,'(A,2ES12.3)') 'exported curve, worst ln c and B errors over the fractions: ', err_c, err_b
            call assert_true(err_c < 0.02_dp, 'the fractions take the contrast of the series')
            call assert_true(err_b < 1.5_dp,  'the fractions take the B-factors of the series')
            call assert_true(abs(log(c(NFRAC-1)/c(NFRAC)) - GAMMA_TRUE) < 0.005_dp, &
                &'past the last window the contrast falls by gamma per unit dose')
            call assert_true(abs((b(NFRAC) - b(NFRAC-1)) - BETA_TRUE) < 0.5_dp, &
                &'past the last window the B-factor rises by beta per unit dose')
        endif
        open(newunit=funit, file=FNAME, status='old')
        close(funit, status='delete')
        call curve%kill
        call readback%kill
        deallocate(em, om, shell, kfreq)
    end subroutine test_exports_the_dose_curve

    !> The per-shell critical exposure assumes no model: on a series faded by Grant and Grigorieff's curve, which
    !! is not one B-factor, it recovers their Ne(k), read back from the decay table the program writes
    subroutine test_shell_exposure()
        integer, parameter :: NW = 5
        character(len=*), parameter :: FNAME = 'dose_series_tester_decay.txt'
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: atten(NW,NSH), sigma(NW), dose(NW), vals(6+4*NW), dev(NSH)
        integer(kind=8) :: seed
        integer         :: w, s, funit, ios, nrow, ngg
        write(*,'(A)') 'test_shell_exposure'
        call shell_grid(shell, kfreq)
        do w = 1, NW
            dose(w) = 4._dp * real(w,dp) - 2._dp
        enddo
        do s = 1, NSH
            do w = 1, NW
                atten(w,s) = gg_atten(dose(w), dose(1), kfreq(s))
            enddo
        enddo
        sigma = 0.01_dp
        seed  = 8080_8
        call synth_atten(NW, atten, sigma, shell, seed, em, om)
        call ds%new(NW, kfreq, 300._dp, .true.)
        do w = 1, NW
            call ds%set_window(w, dose(w) - 2._dp, dose(w) + 2._dp)
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%write_decay_table(string(FNAME))
        call ds%kill
        open(newunit=funit, file=FNAME, status='old', action='read', iostat=ios)
        call assert_true(ios == 0, 'the decay table is written')
        read(funit,*); read(funit,*)
        nrow = 0; ngg = 0
        do
            read(funit,*,iostat=ios) vals
            if( ios /= 0 ) exit
            ! their curve is written as given
            if( abs(vals(6) - (0.245_dp*vals(2)**(-1.665_dp) + 2.81_dp)) < 1.e-6_dp*vals(6) ) ngg = ngg + 1
            if( vals(3) /= vals(3) ) cycle                               ! no measured decay in this shell
            nrow = nrow + 1
            dev(nrow) = abs(vals(3) - vals(6)) / vals(6)
        enddo
        close(funit, status='delete')
        call assert_true(ngg >= 6, 'the Grant and Grigorieff column is their curve')
        call assert_true(nrow >= 6, 'the critical exposure is measured in most shells')
        if( nrow >= 1 ) call assert_true(median(dev(1:nrow)) < 0.1_dp, &
            &'the per-shell critical exposure recovers a curve that is not one B-factor')
        deallocate(em, om, shell, kfreq)
    end subroutine test_shell_exposure

    !> The excess scatter about a weighted line counts the fitted slope in its trace as well as in its degrees of
    !! freedom: an intercept-only trace would give 2.5 for the first case, not 5
    subroutine test_line_tau2()
        real(dp) :: d3(3), b3(3), se3(3), d5(5), b5(5), se5(5), x(5,2), wm(5,5), xtwx(2,2), inv(2,2), pm(5,5), trp, q
        integer  :: i
        write(*,'(A)') 'test_line_tau2'
        ! equal weights, residual [1,-2,1] orthogonal to the line: Q = 6, k - 2 = 1, tr(P) = 1
        d3  = [0._dp, 1._dp, 2._dp]
        se3 = 1._dp
        b3  = 2._dp + 3._dp*d3 + [1._dp, -2._dp, 1._dp]
        call assert_true(abs(line_tau2(b3, d3, se3) - 5._dp) < 1.e-12_dp, 'excess scatter about a line in closed form')
        call assert_true(abs(line_tau2(b3 + 7._dp - 0.5_dp*d3, d3, se3) - 5._dp) < 1.e-12_dp, &
            &'adding a line leaves the excess scatter unchanged')
        call assert_true(line_tau2(b3(1:2), d3(1:2), se3(1:2)) == 0._dp, 'two points carry no scatter about their line')
        ! unequal weights against the matrix form (b'Pb - (k-2)) / tr(P)
        d5  = [1._dp, 3._dp, 4._dp, 7._dp, 9._dp]
        se5 = [0.5_dp, 0.8_dp, 1._dp, 1.3_dp, 2._dp]
        b5  = [4.2_dp, 11.9_dp, 17.5_dp, 27.1_dp, 38.6_dp]
        x(:,1) = 1._dp
        x(:,2) = d5
        wm = 0._dp
        do i = 1, 5
            wm(i,i) = 1._dp / se5(i)**2
        enddo
        xtwx = matmul(transpose(x), matmul(wm, x))
        inv  = reshape([xtwx(2,2), -xtwx(2,1), -xtwx(1,2), xtwx(1,1)], [2,2]) / &
            &(xtwx(1,1)*xtwx(2,2) - xtwx(1,2)*xtwx(2,1))
        pm   = wm - matmul(matmul(wm, x), matmul(inv, matmul(transpose(x), wm)))
        trp  = 0._dp
        do i = 1, 5
            trp = trp + pm(i,i)
        enddo
        q = dot_product(b5, matmul(pm, b5))
        call assert_true(abs(line_tau2(b5, d5, se5) - (q - 3._dp)/trp) < 1.e-9_dp * max(1._dp, abs((q - 3._dp)/trp)), &
            &'unequal weights match the matrix form of the moment estimate')
    end subroutine test_line_tau2

    !> Under bfac_model=free a window carrying no signal must not take the weight: its Guinier fit used to
    !! diverge (B near -2000) and hand it the whole band edge. A window that is genuinely sharper than
    !! window 1 is a good fit and keeps its weight.
    subroutine test_failed_fits_are_rejected()
        integer, parameter :: NREP = 10
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w(:,:)
        real(dp)        :: atten(5,NSH), sigma(5), wdead, wsharp
        integer(kind=8) :: seed
        integer         :: irep, iw, klo, khi, kind, s
        logical         :: l_finite
        write(*,'(A)') 'test_failed_fits_are_rejected'
        call shell_grid(shell, kfreq)
        sigma = 0.15_dp
        do kind = 1, 2
            do s = 1, NSH
                do iw = 1, 4
                    atten(iw,s) = exp(-10._dp*(iw-1) * kfreq(s)**2 / 4._dp)
                enddo
                if( kind == 1 )then
                    atten(5,s) = 0._dp                                    ! no signal at all
                else
                    atten(5,s) = exp(150._dp * kfreq(s)**2 / 4._dp)       ! sharper than window 1
                endif
            enddo
            wdead = 0._dp; wsharp = 1._dp; l_finite = .true.
            do irep = 1, NREP
                seed = 4242_8 * int(irep,8) + int(kind,8)
                call synth_atten(5, atten, sigma, shell, seed, em, om)
                call ds%new(5, kfreq, 300._dp, .false.)
                do iw = 1, 5
                    call ds%set_window(iw, real(iw-1,dp), real(iw,dp))
                enddo
                call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
                call ds%set_band(1, 0)
                call ds%fit
                call ds%get_band(klo, khi)
                call ds%get_weights(w, DOSE_WEIGHTS_MODEL)
                if( .not. all(abs(w) <= 1._dp + 1.e-12_dp) ) l_finite = .false.
                if( kind == 1 )then
                    wdead = max(wdead, maxval(w(5,klo:khi)))
                else
                    wsharp = min(wsharp, w(5,khi))
                endif
                call ds%kill
            enddo
            call assert_true(l_finite, 'model weights stay finite and within [0,1]')
            if( kind == 1 )then
                call assert_true(wdead < 0.05_dp, 'a window without signal takes no weight')
            else
                call assert_true(wsharp > 0.5_dp, 'a window sharper than window 1 keeps its weight at the band edge')
            endif
        enddo
        deallocate(em, om, shell, kfreq, w)
    end subroutine test_failed_fits_are_rejected

    !> The reported error of a_w(s) matches its spread over noise seeds, at high and low SNR. The spread of
    !! the cross products would add the structure's own power fluctuation (errors too large by about the
    !! amplitude SNR); weights alone would miss the noise every pair of a window shares (too small).
    subroutine test_errors_are_calibrated()
        integer,  parameter :: NREP = 24, NW = 4, KB = 20, NSIG = 2
        real(dp), parameter :: SIGS(NSIG) = [0.02_dp, 0.3_dp]
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), a(:,:), se(:,:), grms(:)
        real(dp) :: cw(NW), bw(NW), sigma(NW), sa(NW,KB), sa2(NW,KB), sse(NW,KB), sg(NW), ratio(KB), sd
        integer(kind=8) :: seed
        integer  :: isig, irep, w, s
        write(*,'(A)') 'test_errors_are_calibrated'
        call shell_grid(shell, kfreq)
        cw = 1._dp
        bw = [0._dp, 15._dp, 30._dp, 45._dp]
        do isig = 1, NSIG
            sigma = SIGS(isig)
            sa = 0._dp; sa2 = 0._dp; sse = 0._dp; sg = 0._dp
            do irep = 1, NREP
                seed = 1000_8 * int(isig,8) + int(irep,8)
                call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
                call ds%new(NW, kfreq, 300._dp, .true.)
                do w = 1, NW
                    call ds%set_window(w, real(w-1,dp), real(w,dp))
                enddo
                call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
                call ds%set_band(1, KB)
                call ds%fit
                call ds%get_attenuation(a, se, grms)
                sa  = sa  + a(:,1:KB)
                sa2 = sa2 + a(:,1:KB)**2
                sse = sse + se(:,1:KB)
                sg  = sg  + grms
                call ds%kill
            enddo
            do w = 2, NW
                do s = 1, KB
                    sd       = sqrt(max(1.e-30_dp, (sa2(w,s) - sa(w,s)**2/NREP) / (NREP - 1)))
                    ratio(s) = (sse(w,s)/NREP) / sd
                enddo
                call assert_true(median(ratio) > 0.8_dp .and. median(ratio) < 1.25_dp, &
                    &'reported attenuation error matches its spread over seeds')
                ! a two-parameter Guinier fit over KB shells leaves sqrt((KB-2)/KB) of rms
                call assert_true(sg(w)/NREP > 0.75_dp .and. sg(w)/NREP < 1.15_dp, &
                    &'the Guinier rms of a series on the model is about one')
            enddo
        enddo
        deallocate(em, om, shell, kfreq, a, se, grms)
    end subroutine test_errors_are_calibrated

    !> the independent fraction of a shell under a mask: one without a mask, and 3/(1 + 2N/t) for a slab of
    !! thickness t, whose projections are a constant along its normal and a slab across it
    subroutine test_shell_independence()
        integer, parameter :: N = 16, T = 4
        real    :: m(N,N,N)
        real(dp) :: ball, fvol
        integer :: i, j, k
        write(*,'(A)') 'test_shell_independence'
        m = 1.
        call assert_true(abs(shell_independence(m) - 1._dp) < 1.e-12_dp, 'no mask leaves every coefficient independent')
        m = 0.
        m(:,:,5:4+T) = 1.
        call assert_true(abs(shell_independence(m) - 3._dp/(1._dp + 2._dp*real(N,dp)/real(T,dp))) < 1.e-12_dp, &
            &'a slab leaves the closed-form fraction')
        ! a ball follows the plane-section law, about 1.075 times its volume fraction to the 2/3
        m = 0.
        do k = 1, N
            do j = 1, N
                do i = 1, N
                    if( (i-N/2-1)**2 + (j-N/2-1)**2 + (k-N/2-1)**2 <= 25 ) m(i,j,k) = 1.
                enddo
            enddo
        enddo
        fvol = real(count(m > 0.5), dp) / real(N**3, dp)
        ball = shell_independence(m)
        call assert_true(abs(ball / (1.075_dp * fvol**(2._dp/3._dp)) - 1._dp) < 0.1_dp, &
            &'a ball follows the plane-section law')
        call assert_true(ball > fvol, 'a shell keeps more independence than the volume fraction suggests')
        ! across shells: the whole kernel against its plane, (4/3) N/D for a ball of diameter D
        call assert_true(abs(shell_correlation(m) / (4._dp*real(N,dp)/(3._dp*10._dp)) - 1._dp) < 0.15_dp, &
            &'a ball inflates a fit across shells by about (4/3) N/D')
        m = 1.
        call assert_true(abs(shell_correlation(m) - 1._dp) < 1.e-12_dp, 'no mask leaves the shells uncorrelated')
        m = 0.
        m(:,:,5:4+T) = 1.
        call assert_true(abs(shell_correlation(m) - (real(N,dp)/real(T,dp) + 2._dp)/3._dp) < 1.e-12_dp, &
            &'a slab correlates shells along its thin axis alone')
    end subroutine test_shell_independence

    !> A weak reference correlates every relative attenuation: exposure errors must match independent seeds.
    subroutine test_exposure_errors_share_the_reference()
        use simple_syslib, only: del_file
        integer, parameter :: NW = 10, NREP = 64, KSHELL = 18
        character(len=*), parameter :: FNAME = 'dose_series_tester_exposure_errors.txt'
        character(len=*), parameter :: STEPS = 'dose_series_tester_step_errors.txt'
        type(dose_series) :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer, allocatable :: shell(:)
        real(dp), allocatable :: kfreq(:)
        real(dp) :: cw(NW), bw(NW), sigma(NW), row(7), values(NREP), errors(NREP), step_values(NREP), step_errors(NREP), truth, ratio
        integer(kind=8) :: seed
        integer :: rep, w, u, s
        write(*,'(A)') 'test_exposure_errors_share_the_reference'
        call shell_grid(shell, kfreq)
        cw = 1._dp
        bw = [(10._dp*real(w-1,dp), w=1,NW)]
        sigma = 0.01_dp
        sigma(1) = 0.08_dp
        truth = 2._dp / (10._dp*kfreq(KSHELL)**2)
        do rep = 1, NREP
            seed = 100000_8 + int(rep,8)*101_8
            call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
            call ds%new(NW, kfreq, 300._dp, .true.)
            do w = 1, NW
                call ds%set_window(w, real(w-1,dp), real(w,dp))
            enddo
            call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
            call ds%set_band(1, NSH)
            call ds%fit
            call ds%write_decay_table(string(FNAME))
            open(newunit=u, file=FNAME, status='old', action='read')
            read(u,*); read(u,*)
            do s = 1, KSHELL
                read(u,*) row
            enddo
            close(u)
            values(rep) = row(3)
            errors(rep) = row(4)
            call ds%write_step_exposure(string(STEPS), [1._dp/kfreq(KSHELL)])
            open(newunit=u, file=STEPS, status='old', action='read')
            read(u,*)
            do w = 1, 5
                read(u,*) row(:5)
            enddo
            close(u)
            step_values(rep) = row(3)
            step_errors(rep) = row(4)
            call ds%kill
        enddo
        ratio = sqrt(sum((values-sum(values)/NREP)**2)/real(NREP-1,dp)) / (sum(errors)/NREP)
        write(*,'(A,F9.4)') 'critical-exposure scatter / mean reported error: ', ratio
        call assert_true(abs(sum(values)/NREP-truth) < 0.02_dp*truth, 'critical exposure recovers its analytic truth')
        call assert_true(ratio > 0.75_dp .and. ratio < 1.35_dp, 'exposure errors include the common noisy reference')
        ratio = sqrt(sum((step_values-sum(step_values)/NREP)**2)/real(NREP-1,dp)) / (sum(step_errors)/NREP)
        write(*,'(A,F9.4)') 'neighboring-window exposure scatter / mean reported error: ', ratio
        call assert_true(abs(sum(step_values)/NREP-truth) < 0.02_dp*truth, 'step exposure recovers its analytic truth')
        call assert_true(ratio > 0.75_dp .and. ratio < 1.35_dp, 'step exposure errors cancel common-reference noise')
        call del_file(FNAME)
        call del_file(STEPS)
    end subroutine test_exposure_errors_share_the_reference

    !> Differences of free B-factors cancel the common reference noise; compare with independent realizations.
    subroutine test_relative_b_errors_share_the_reference()
        use simple_syslib, only: del_file
        integer, parameter :: NW = 5, NREP = 48
        character(len=*), parameter :: FNAME = 'dose_series_tester_b_errors.txt'
        type(dose_series) :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer, allocatable :: shell(:)
        real(dp), allocatable :: kfreq(:)
        real(dp) :: cw(NW), bw(NW), sigma(NW), row(10), values(NREP), errors(NREP), ratio
        integer(kind=8) :: seed
        integer :: rep, w, u
        write(*,'(A)') 'test_relative_b_errors_share_the_reference'
        call shell_grid(shell, kfreq)
        cw = 1._dp
        bw = [0._dp, -6._dp, 10._dp, 30._dp, 50._dp]
        sigma = 0.01_dp
        sigma(1) = 0.08_dp
        do rep = 1, NREP
            seed = 200000_8 + int(rep,8)*101_8
            call synth(NW, cw, bw, sigma, kfreq, shell, seed, em, om)
            call ds%new(NW, kfreq, 300._dp, .false.)
            do w = 1, NW
                call ds%set_window(w, real(w-1,dp), real(w,dp))
            enddo
            call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
            call ds%set_band(1, NSH)
            call ds%fit
            call ds%write_model_table(string(FNAME))
            open(newunit=u, file=FNAME, status='old', action='read')
            read(u,*); read(u,*)
            do w = 1, 3
                read(u,*) row
            enddo
            close(u)
            values(rep) = row(4)
            errors(rep) = row(5)
            call ds%kill
        enddo
        ratio = sqrt(sum((values-sum(values)/NREP)**2)/real(NREP-1,dp)) / (sum(errors)/NREP)
        write(*,'(A,F9.4)') 'relative B scatter / mean reported error: ', ratio
        call assert_true(abs(sum(values)/NREP-16._dp) < 0.1_dp, 'relative B recovers the planted difference')
        call assert_true(ratio > 0.75_dp .and. ratio < 1.35_dp, 'relative B error includes the shared-reference covariance')
        call del_file(FNAME)
    end subroutine test_relative_b_errors_share_the_reference

    !> Comparison tables expose paired error differences without an uncalibrated significance claim.
    subroutine check_descriptive_cv_table( ds )
        use simple_syslib, only: del_file
        class(dose_series), intent(in) :: ds
        character(len=*), parameter :: FNAME = 'dose_series_tester_cv.txt'
        character(len=512) :: text
        integer :: u
        call ds%write_cv_table(string(FNAME))
        open(newunit=u, file=FNAME, status='old', action='read')
        read(u,'(A)') text
        call assert_true(index(text, 'descriptive') > 0, 'CV header identifies descriptive comparisons')
        read(u,'(A)') text
        call assert_true(index(text, 'error_model_minus_control') > 0 .and. index(text, 'z_model') == 0, &
            &'the paired output is an error difference, not a z-score')
        close(u)
        call del_file(FNAME)
    end subroutine check_descriptive_cv_table

    real(dp) function median( x ) result( m )
        real(dp), intent(in) :: x(:)
        real(dp) :: y(size(x)), t
        integer  :: i, j
        y = x
        do i = 2, size(y)
            t = y(i)
            j = i - 1
            do while( j >= 1 )
                if( y(j) <= t ) exit
                y(j+1) = y(j)
                j = j - 1
            enddo
            y(j+1) = t
        enddo
        m = y((size(y)+1)/2)
    end function median

    !> a deterministic generator, so the suite never depends on global RNG state
    real(dp) function next_u( state ) result( u )
        integer(kind=8), intent(inout) :: state
        state = modulo(1103515245_8 * state + 12345_8, 2147483648_8)
        u     = real(state, dp) / 2147483648._dp
    end function next_u

    !> one standard complex Gaussian, unit total power
    complex(dp) function next_z( state ) result( z )
        integer(kind=8), intent(inout) :: state
        real(dp) :: u1, u2, r
        u1 = max(1.e-12_dp, next_u(state))
        u2 = next_u(state)
        r  = sqrt(-log(u1))
        z  = cmplx(r*cos(2._dp*PI*u2), r*sin(2._dp*PI*u2), kind=dp)
    end function next_z

    subroutine shell_grid( shell, kfreq )
        integer,  allocatable, intent(inout) :: shell(:)
        real(dp), allocatable, intent(inout) :: kfreq(:)
        integer :: s, i, j
        if( allocated(shell) ) deallocate(shell)
        if( allocated(kfreq) ) deallocate(kfreq)
        allocate(shell(NSH*NPERSH), kfreq(NSH))
        i = 0
        do s = 1, NSH
            do j = 1, NPERSH
                i = i + 1
                shell(i) = s
            enddo
            kfreq(s) = real(s, dp) / (real(BOX, dp) * SMPD)
        enddo
    end subroutine shell_grid

    !> Grant & Grigorieff attenuation of a window at dose d against one at dose d1, at frequency k
    real(dp) function gg_atten( d, d1, k ) result( a )
        real(dp), intent(in) :: d, d1, k
        a = exp(-(d - d1) / (2._dp * (0.245_dp * k**(-1.665_dp) + 2.81_dp)))
    end function gg_atten

    !> generate a series under the B-factor model, with independent half noise
    subroutine synth( nwin, cw, bw, sigma, kfreq, shell, seed, em, om )
        integer,                  intent(in)    :: nwin
        real(dp),                 intent(in)    :: cw(:), bw(:), sigma(:), kfreq(:)
        integer,                  intent(in)    :: shell(:)
        integer(kind=8),          intent(inout) :: seed
        complex(sp), allocatable, intent(inout) :: em(:,:), om(:,:)
        real(dp) :: atten(nwin,size(kfreq))
        integer  :: w, s
        do s = 1, size(kfreq)
            do w = 1, nwin
                atten(w,s) = cw(w) * exp(-bw(w) * kfreq(s)**2 / 4._dp)
            enddo
        enddo
        call synth_atten(nwin, atten, sigma, shell, seed, em, om)
    end subroutine synth

    !> generate a series under any attenuation, with independent half noise
    subroutine synth_atten( nwin, atten, sigma, shell, seed, em, om )
        integer,                  intent(in)    :: nwin
        real(dp),                 intent(in)    :: atten(:,:), sigma(:)
        integer,                  intent(in)    :: shell(:)
        integer(kind=8),          intent(inout) :: seed
        complex(sp), allocatable, intent(inout) :: em(:,:), om(:,:)
        complex(dp) :: sig
        real(dp)    :: pwr
        integer     :: i, w, s, ncoeff
        ncoeff = size(shell)
        if( allocated(em) ) deallocate(em)
        if( allocated(om) ) deallocate(om)
        allocate(em(ncoeff,nwin), om(ncoeff,nwin))
        do i = 1, ncoeff
            s   = shell(i)
            pwr = exp(-real(s,dp)/10._dp)            ! protein-like falling shell power
            sig = sqrt(pwr) * next_z(seed)
            do w = 1, nwin
                em(i,w) = cmplx(atten(w,s)*sig + sigma(w)*next_z(seed), kind=sp)
                om(i,w) = cmplx(atten(w,s)*sig + sigma(w)*next_z(seed), kind=sp)
            enddo
        enddo
    end subroutine synth_atten

    !> with equal, vanishing noise the weights follow the attenuation alone; exact zero noise is
    !! read as a dead shell, the safe convention on real data, so the limit is taken from above
    subroutine test_noise_free_attenuation()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w(:,:)
        real(dp)        :: cw(3), bw(3), sigma(3), truth, got
        integer(kind=8) :: seed
        integer         :: s, klo, khi
        write(*,'(A)') 'test_noise_free_attenuation'
        call shell_grid(shell, kfreq)
        cw    = [1._dp, 1._dp, 1._dp]
        bw    = [0._dp, 30._dp, 60._dp]
        sigma = [1.e-3_dp, 1.e-3_dp, 1.e-3_dp]
        seed  = 20261104_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 1._dp)
        call ds%set_window(2, 1._dp, 2._dp)
        call ds%set_window(3, 2._dp, 3._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%get_band(klo, khi)
        call assert_true(khi >= 6, 'noise-free band spans at least six shells')
        call ds%get_weights(w, DOSE_WEIGHTS_MODEL)
        ! the weights of a noiseless series follow the attenuation alone
        do s = klo, khi
            truth = exp(-bw(3) * kfreq(s)**2 / 4._dp) / &
                &(exp(-bw(1)*kfreq(s)**2/4._dp) + exp(-bw(2)*kfreq(s)**2/4._dp) + exp(-bw(3)*kfreq(s)**2/4._dp))
            got   = w(3,s)
            call assert_true(abs(got - truth) < TOL_LOOSE, 'noise-free weight follows the attenuation')
        enddo
        call ds%kill
        deallocate(em, om, shell, kfreq, w)
    end subroutine test_noise_free_attenuation

    !> beta is recovered from a series built on a line in dose
    subroutine test_recovers_the_dose_curve()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(5), bw(5), sigma(5), dose(5), beta_true
        integer(kind=8) :: seed
        integer         :: w
        write(*,'(A)') 'test_recovers_the_dose_curve'
        call shell_grid(shell, kfreq)
        beta_true = 4._dp
        do w = 1, 5
            dose(w)  = real(w, dp) - 0.5_dp
            bw(w)    = beta_true * (dose(w) - dose(1))
            cw(w)    = 1._dp
            sigma(w) = 0.02_dp
        enddo
        seed = 777_8
        call synth(5, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(5, kfreq, 300._dp, .true.)
        do w = 1, 5
            call ds%set_window(w, dose(w) - 0.5_dp, dose(w) + 0.5_dp)
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call assert_true(abs(ds%get_beta() - beta_true) < 1._dp, 'beta recovered from a series on a line')
        call ds%kill
        deallocate(em, om, shell, kfreq)
    end subroutine test_recovers_the_dose_curve

    !> with no damage the matched filter is the precision-weighted mean
    subroutine test_no_damage_gives_precision_weights()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w(:,:)
        real(dp)        :: cw(3), bw(3), sigma(3), expect, den
        integer(kind=8) :: seed
        integer         :: s, klo, khi, v
        write(*,'(A)') 'test_no_damage_gives_precision_weights'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = 0._dp
        sigma = [0.05_dp, 0.10_dp, 0.20_dp]
        seed  = 31337_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 1._dp)
        call ds%set_window(2, 1._dp, 2._dp)
        call ds%set_window(3, 2._dp, 3._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%get_band(klo, khi)
        call ds%get_weights(w, DOSE_WEIGHTS_MODEL)
        do s = klo, khi
            den = 0._dp
            do v = 1, 3
                den = den + 1._dp/sigma(v)**2
            enddo
            expect = (1._dp/sigma(1)**2) / den
            call assert_true(abs(w(1,s) - expect) < TOL_LOOSE, 'undamaged weights follow precision alone')
        enddo
        call ds%kill
        deallocate(em, om, shell, kfreq, w)
    end subroutine test_no_damage_gives_precision_weights

    !> two windows are the smallest identifiable series
    subroutine test_two_windows_are_identifiable()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(2), bw(2), sigma(2)
        integer(kind=8) :: seed
        write(*,'(A)') 'test_two_windows_are_identifiable'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 20._dp]
        sigma = [0.02_dp, 0.02_dp]
        seed  = 99_8
        call synth(2, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(2, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 2._dp)
        call ds%set_window(2, 2._dp, 4._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        ! B rises by 20 A^2 over a dose step of 2, so the line through the gauge has slope 10
        call assert_true(abs(ds%get_beta() - 10._dp) < 4._dp, 'two windows determine the line')
        call ds%kill
        deallocate(em, om, shell, kfreq)
    end subroutine test_two_windows_are_identifiable

    !> the accumulation is exactly the weighted mean of the inputs, coefficient by coefficient
    subroutine test_accumulation_is_a_weighted_mean()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:), ae(:), ao(:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w(:,:)
        real(dp)        :: cw(3), bw(3), sigma(3)
        complex(dp)     :: expect
        integer(kind=8) :: seed
        integer         :: i, v
        write(*,'(A)') 'test_accumulation_is_a_weighted_mean'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 15._dp, 30._dp]
        sigma = [0.05_dp, 0.05_dp, 0.05_dp]
        seed  = 4242_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 1._dp)
        call ds%set_window(2, 1._dp, 2._dp)
        call ds%set_window(3, 2._dp, 3._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%get_weights(w, DOSE_WEIGHTS_MODEL)
        allocate(ae(size(shell)), ao(size(shell)))
        call ds%accumulate(em, om, shell, DOSE_WEIGHTS_MODEL, ae, ao)
        do i = 1, size(shell), 97
            expect = (0._dp, 0._dp)
            do v = 1, 3
                expect = expect + w(v,shell(i)) * cmplx(em(i,v), kind=dp)
            enddo
            call assert_true(abs(cmplx(ae(i), kind=dp) - expect) < 1.e-5_dp, &
                &'accumulated even coefficient is the weighted mean')
        enddo
        call assert_true(sum(abs(w(:,1))) > 0._dp, 'weights are populated outside the band too')
        call ds%kill
        deallocate(em, om, ae, ao, shell, kfreq, w)
    end subroutine test_accumulation_is_a_weighted_mean

    !> the band reaches as far as the best window, not as far as the flat sum
    subroutine test_band_follows_the_best_window()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(3), bw(3), sigma(3)
        integer(kind=8) :: seed
        integer         :: klo, khi, klo2, khi2
        write(*,'(A)') 'test_band_follows_the_best_window'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 200._dp, 400._dp]      ! windows 2 and 3 are dead at high resolution
        sigma = [0.02_dp, 0.02_dp, 0.02_dp]
        seed  = 5150_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 1._dp)
        call ds%set_window(2, 1._dp, 2._dp)
        call ds%set_window(3, 2._dp, 3._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%get_band(klo, khi)
        call assert_int(1, klo, 'band starts at the first shell')
        call assert_true(khi > NSH/2, 'band follows the least damaged window')
        ! an explicit lp overrides the derived band
        call ds%set_band(1, 10)
        call ds%get_band(klo2, khi2)
        call assert_int(10, khi2, 'an explicit lp shell sets the band')
        ! a high-pass shell moves the start of the band and nothing else
        call ds%set_band(4, 0)
        call ds%get_band(klo2, khi2)
        call assert_int(4, klo2, 'a high-pass shell sets the start of the band')
        call assert_int(khi, khi2, 'the high-pass leaves the derived end of the band alone')
        call ds%kill
        deallocate(em, om, shell, kfreq)
    end subroutine test_band_follows_the_best_window

    !> a global scale on the coefficients, such as an FFT normalisation, changes nothing
    subroutine test_statistics_are_scale_free()
        type(dose_series)        :: ds1, ds2
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w1(:,:), w2(:,:)
        real(dp)        :: cw(4), bw(4), sigma(4)
        real(sp), parameter :: SCALE = 1.e-4
        integer(kind=8) :: seed
        integer         :: klo1, khi1, klo2, khi2, w
        write(*,'(A)') 'test_statistics_are_scale_free'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 10._dp, 20._dp, 30._dp]
        sigma = [0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp]
        seed  = 1234_8
        call synth(4, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds1%new(4, kfreq, 300._dp, .true.)
        call ds2%new(4, kfreq, 300._dp, .true.)
        do w = 1, 4
            call ds1%set_window(w, real(w-1,dp), real(w,dp))
            call ds2%set_window(w, real(w-1,dp), real(w,dp))
        enddo
        call ds1%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds2%accumulate_gram(SCALE*em, SCALE*om, shell, 1._dp, 1._dp)
        call ds1%set_band(1, 0)
        call ds2%set_band(1, 0)
        call ds1%get_band(klo1, khi1)
        call ds2%get_band(klo2, khi2)
        call assert_int(khi1, khi2, 'the derived band does not depend on the coefficient scale')
        call ds1%fit
        call ds2%fit
        call assert_true(abs(ds1%get_beta() - ds2%get_beta()) < 1.e-6_dp*max(1._dp, abs(ds1%get_beta())), &
            &'beta does not depend on the coefficient scale')
        call ds1%get_weights(w1, DOSE_WEIGHTS_MODEL)
        call ds2%get_weights(w2, DOSE_WEIGHTS_MODEL)
        call assert_true(maxval(abs(w1 - w2)) < 1.e-6_dp, 'the weights do not depend on the coefficient scale')
        call ds1%kill
        call ds2%kill
        deallocate(em, om, shell, kfreq, w1, w2)
    end subroutine test_statistics_are_scale_free

    !> the control's weights are the matched filter under Grant & Grigorieff and the flat mean's under
    !! a = 1, so their ratio is the GG attenuation whatever precisions the data gave
    subroutine test_control_and_flat_weights()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), wc(:,:), wf(:,:)
        real(dp)        :: cw(3), bw(3), sigma(3), dose(3), ratio, truth
        integer(kind=8) :: seed
        integer         :: s, w, klo, khi
        write(*,'(A)') 'test_control_and_flat_weights'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 15._dp, 30._dp]
        sigma = [0.05_dp, 0.08_dp, 0.12_dp]
        dose  = [0.5_dp, 1.5_dp, 2.5_dp]
        seed  = 2718_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        do w = 1, 3
            call ds%set_window(w, dose(w) - 0.5_dp, dose(w) + 0.5_dp)
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%get_band(klo, khi)
        call ds%get_weights(wc, DOSE_WEIGHTS_CONTROL)
        call ds%get_weights(wf, DOSE_WEIGHTS_FLAT)
        do s = klo, khi
            call assert_true(abs(sum(wc(:,s)) - 1._dp) < 1.e-12_dp, 'control weights sum to one')
            call assert_true(abs(sum(wf(:,s)) - 1._dp) < 1.e-12_dp, 'flat weights sum to one')
            do w = 2, 3
                ratio = (wc(w,s)/wf(w,s)) / (wc(1,s)/wf(1,s))
                truth = gg_atten(dose(w), dose(1), kfreq(s))
                call assert_true(abs(ratio - truth) < 1.e-10_dp, 'control over flat weight is the GG attenuation')
            enddo
        enddo
        call ds%kill
        deallocate(em, om, shell, kfreq, wc, wf)
    end subroutine test_control_and_flat_weights

    !> the zero-dose map is one real factor per shell, the band edge's factor held past the band,
    !! and about one where nothing was damaged
    subroutine test_zerodose_rescaling()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:), s0(:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real(dp)        :: cw(3), bw(3), sigma(3), fac(NSH), r
        logical         :: l_const, l_held, l_near_one
        integer(kind=8) :: seed
        integer         :: i, s, klo, khi
        write(*,'(A)') 'test_zerodose_rescaling'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = 0._dp
        sigma = [0.03_dp, 0.03_dp, 0.03_dp]
        seed  = 1618_8
        call synth(3, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(3, kfreq, 300._dp, .true.)
        call ds%set_window(1, 0._dp, 1._dp)
        call ds%set_window(2, 1._dp, 2._dp)
        call ds%set_window(3, 2._dp, 3._dp)
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 10)
        call ds%fit
        call ds%get_band(klo, khi)
        allocate(s0(size(shell)))
        call ds%zerodose(em(:,1), shell, s0)
        fac     = 0._dp
        l_const = .true.
        do i = 1, size(shell)
            if( abs(em(i,1)) < 1.e-6 ) cycle
            r = real(cmplx(s0(i), kind=dp) / cmplx(em(i,1), kind=dp), dp)
            s = shell(i)
            if( fac(s) == 0._dp )then
                fac(s) = r
            else if( abs(r - fac(s)) > 1.e-4_dp * abs(fac(s)) )then
                l_const = .false.
            endif
        enddo
        call assert_true(l_const, 'the zero-dose map is one factor per shell')
        l_held = .true.
        do s = khi + 1, NSH
            if( abs(fac(s) - fac(khi)) > 1.e-4_dp * abs(fac(khi)) ) l_held = .false.
        enddo
        call assert_true(l_held, 'past the band the band edge factor is held')
        l_near_one = .true.
        do s = klo, khi
            if( abs(fac(s) - 1._dp) > TOL_LOOSE ) l_near_one = .false.
        enddo
        call assert_true(l_near_one, 'with no damage the zero-dose factor is about one')
        call ds%kill
        deallocate(em, om, s0, shell, kfreq)
    end subroutine test_zerodose_rescaling

    !> an undamaged series makes the Grant & Grigorieff control mispredict its late windows, so the
    !! paired mean difference must be negative where the control's attenuation bites
    subroutine test_paired_cv_rejects_a_wrong_control()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), cv(:,:), delta(:)
        real(dp)        :: cw(5), bw(5), sigma(5), cvm, cvc
        integer(kind=8) :: seed
        integer         :: w, s, klo, khi, nwin_model, nupper
        write(*,'(A)') 'test_paired_cv_rejects_a_wrong_control'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = 0._dp
        sigma = 0.03_dp
        seed  = 8128_8
        call synth(5, cw, bw, sigma, kfreq, shell, seed, em, om)
        call ds%new(5, kfreq, 300._dp, .true.)
        do w = 1, 5
            call ds%set_window(w, real(w-1,dp), real(w,dp))
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%paired_cv(em, om, shell)
        call ds%get_band(klo, khi)
        call ds%get_cv(cv, delta)
        call assert_true(maxval(abs(delta(klo:khi) - &
            &(cv(DOSE_WEIGHTS_MODEL,klo:khi)-cv(DOSE_WEIGHTS_CONTROL,klo:khi)))) < 1.e-9_dp, &
            &'coefficient-paired differences equal independently computed Gram error differences')
        call check_descriptive_cv_table(ds)
        cvm = 0._dp; cvc = 0._dp; nwin_model = 0; nupper = 0
        do s = (klo + khi)/2, khi
            cvm    = cvm + cv(DOSE_WEIGHTS_MODEL,s)
            cvc    = cvc + cv(DOSE_WEIGHTS_CONTROL,s)
            nupper = nupper + 1
            if( delta(s) < 0._dp ) nwin_model = nwin_model + 1
        enddo
        call assert_true(cvm < cvc, 'the model predicts an undamaged series better than the control')
        call assert_true(2*nwin_model >= nupper, 'the paired delta favours the model in most upper shells')
        call ds%kill
        deallocate(em, om, shell, kfreq, cv, delta)
    end subroutine test_paired_cv_rejects_a_wrong_control

    !> damage of the Grant & Grigorieff shape cannot be predicted better by one B-factor per window
    subroutine test_paired_cv_keeps_a_true_control()
        type(dose_series)        :: ds
        complex(sp), allocatable :: em(:,:), om(:,:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), cv(:,:), delta(:)
        real(dp)        :: atten(5,NSH), sigma(5), dose(5), cvm, cvc
        integer(kind=8) :: seed
        integer         :: w, s, klo, khi
        write(*,'(A)') 'test_paired_cv_keeps_a_true_control'
        call shell_grid(shell, kfreq)
        do w = 1, 5
            dose(w) = real(w, dp) - 0.5_dp
        enddo
        do s = 1, NSH
            do w = 1, 5
                atten(w,s) = gg_atten(dose(w), dose(1), kfreq(s))
            enddo
        enddo
        sigma = 0.01_dp
        seed  = 496_8
        call synth_atten(5, atten, sigma, shell, seed, em, om)
        call ds%new(5, kfreq, 300._dp, .true.)
        do w = 1, 5
            call ds%set_window(w, dose(w) - 0.5_dp, dose(w) + 0.5_dp)
        enddo
        call ds%accumulate_gram(em, om, shell, 1._dp, 1._dp)
        call ds%set_band(1, 0)
        call ds%fit
        call ds%paired_cv(em, om, shell)
        call ds%get_band(klo, khi)
        call ds%get_cv(cv, delta)
        cvm = sum(cv(DOSE_WEIGHTS_MODEL,klo:khi))
        cvc = sum(cv(DOSE_WEIGHTS_CONTROL,klo:khi))
        call assert_true(cvc <= cvm, 'the true curve predicts at least as well as the fitted model')
        call ds%kill
        deallocate(em, om, shell, kfreq, cv, delta)
    end subroutine test_paired_cv_keeps_a_true_control

    !> every parallel loop sums each shell in one fixed order, so a team reproduces one thread exactly
    subroutine test_parallel_matches_serial()
        complex(sp), allocatable :: em(:,:), om(:,:), ae1(:), ao1(:), ae2(:), ao2(:), s01(:), s02(:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:), w1(:,:), w2(:,:), cv1(:,:), cv2(:,:), delta1(:), delta2(:)
        real(dp)        :: cw(4), bw(4), sigma(4)
        integer(kind=8) :: seed
        integer         :: nthr_prev, nteam
        write(*,'(A)') 'test_parallel_matches_serial'
        call shell_grid(shell, kfreq)
        cw    = 1._dp
        bw    = [0._dp, 12._dp, 24._dp, 36._dp]
        sigma = 0.05_dp
        seed  = 6174_8
        call synth(4, cw, bw, sigma, kfreq, shell, seed, em, om)
        nthr_prev = 1
        nteam     = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ nteam     = max(2, min(4, omp_get_num_procs()))
        !$ call omp_set_num_threads(1)
        call run_pipeline(em, om, shell, kfreq, w1, cv1, delta1, ae1, ao1, s01)
        !$ call omp_set_num_threads(nteam)
        call run_pipeline(em, om, shell, kfreq, w2, cv2, delta2, ae2, ao2, s02)
        !$ call omp_set_num_threads(nthr_prev)
        call assert_true(all(w1 == w2),   'weights at a team equal those at one thread')
        call assert_true(all(cv1 == cv2), 'held-out errors at a team equal those at one thread')
        call assert_true(all(delta1 == delta2),   'paired delta at a team equals that at one thread')
        call assert_true(all(ae1 == ae2) .and. all(ao1 == ao2), 'accumulations at a team equal those at one thread')
        call assert_true(all(s01 == s02), 'the zero-dose map at a team equals that at one thread')
        deallocate(em, om, shell, kfreq, w1, w2, cv1, cv2, delta1, delta2, ae1, ao1, ae2, ao2, s01, s02)
    end subroutine test_parallel_matches_serial

    subroutine run_pipeline( em, om, shell, kfreq, w, cv, delta, ae, ao, s0 )
        complex(sp),              intent(in)    :: em(:,:), om(:,:)
        integer,                  intent(in)    :: shell(:)
        real(dp),                 intent(in)    :: kfreq(:)
        real(dp),    allocatable, intent(inout) :: w(:,:), cv(:,:), delta(:)
        complex(sp), allocatable, intent(inout) :: ae(:), ao(:), s0(:)
        type(dose_series) :: ds
        integer :: iw
        call ds%new(size(em,2), kfreq, 300._dp, .true.)
        do iw = 1, size(em,2)
            call ds%set_window(iw, real(iw-1,dp), real(iw,dp))
        enddo
        call ds%accumulate_gram(em, om, shell, 0.5_dp, 1._dp)
        call ds%set_band(2, 0)
        call ds%fit
        call ds%paired_cv(em, om, shell)
        call ds%get_weights(w, DOSE_WEIGHTS_MODEL)
        call ds%get_cv(cv, delta)
        allocate(ae(size(shell)), ao(size(shell)), s0(size(shell)))
        call ds%accumulate(em, om, shell, DOSE_WEIGHTS_MODEL, ae, ao)
        call ds%zerodose(ae, shell, s0)
        call ds%kill
    end subroutine run_pipeline

end module simple_dose_series_tester
