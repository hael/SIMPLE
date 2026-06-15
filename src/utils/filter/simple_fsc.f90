!@descr: various Fourier Shell Correlation utilities
module simple_fsc
use simple_core_module_api
use simple_image, only: image
use CPlot2D_wrapper_module
implicit none

public :: phase_rand_fsc, plot_fsc, plot_fsc2, plot_phrand_fsc, plot_fsc_area_score
public :: fsc_area_score_result, write_fsc_area_score_outputs
private
#include "simple_local_flags.inc"

type :: fsc_area_score_result
    integer :: nfreqs=0
    integer :: ndirs=0
    integer :: min_count=1
    integer :: min_dir=0
    integer :: max_dir=0
    real    :: cone_half_angle_deg=20.
    real    :: threshold=0.143
    real    :: cfar=0.
    real    :: wauc_min=0.
    real    :: wauc_max=0.
    real,    allocatable :: res(:)
    real,    allocatable :: dirs(:,:)
    real,    allocatable :: cfsc(:,:)
    real,    allocatable :: wauc(:)
    real,    allocatable :: crossing_find(:)
    real,    allocatable :: crossing_res(:)
    integer, allocatable :: counts(:,:)
    integer, allocatable :: included_bins(:)
    logical              :: exists = .false.
  contains
    procedure :: new
    procedure :: calc_fsc_area_score
    procedure :: kill
end type fsc_area_score_result

contains

    ! calculate phase-randomized FSC according to Chen et al,JSB,2013
    subroutine phase_rand_fsc(even, odd, envmask, msk, state, n, fsc, fsc_t, fsc_n)
        class(image),            intent(inout) :: even, odd, envmask
        integer,                 intent(in)    :: state, n
        real,                    intent(in)    :: msk
        real, allocatable,       intent(out)   :: fsc(:), fsc_t(:), fsc_n(:)
        real    :: lp_rand, smpd
        integer :: ldim(3), k, k_rand
        call random_seed()
        ldim = even%get_ldim()
        smpd = even%get_smpd()
        allocate(fsc(n),fsc_t(n),fsc_n(n), source=0.)
        call even%ifft()                        ! Fourier space
        call odd%ifft()
        ! Enveloppe-masked FSC
        call even%zero_env_background(envmask)
        call odd%zero_env_background(envmask)
        call even%mul(envmask)                  ! mask
        call odd%mul(envmask)
        call even%fft()                         ! Fourier space
        call odd%fft()
        call even%fsc(odd, fsc_t)               ! FSC
        ! Randomize then calculate masked FSC
        k_rand = get_find_at_crit(fsc_t, ENVMSK_FSC_THRESH)
        if( k_rand > n-3 )then
            ! reverts to sherical masking
            call even%ifft()
            call odd%ifft()
            call even%mask3D_soft(msk)
            call odd%mask3D_soft(msk)
            call even%fft()
            call odd%fft()
            call even%fsc(odd, fsc)
        else
            lp_rand = calc_lowpass_lim(k_rand, ldim(1), smpd)
            ! randomize
            call even%phase_rand(lp_rand)
            call odd%phase_rand(lp_rand)
            ! mask
            call even%ifft()
            call odd%ifft()
            call even%zero_env_background(envmask)
            call odd%zero_env_background(envmask)
            call even%mul(envmask)
            call odd%mul(envmask)
            ! FSC phase-randomized
            call even%fft()
            call odd%fft()
            call even%fsc(odd, fsc_n)
            ! correction
            fsc = fsc_t
            do k = k_rand+2,n
                fsc(k) = (fsc_t(k)-fsc_n(k)) / (1.-fsc_n(k))
            enddo
            call arr2file(fsc_t, string('fsct_state'//int2str_pad(state,2)//BIN_EXT))
            call arr2file(fsc_n, string('fscn_state'//int2str_pad(state,2)//BIN_EXT))
        endif
    end subroutine phase_rand_fsc

    subroutine plot_fsc( n, fsc, res, smpd, tmpl_fname )
        integer,           intent(in) :: n
        real,              intent(in) :: fsc(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        type(string)              :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, fname_eps
        integer  :: k,iostat
        if( n == 0 ) THROW_HARD('Empty FSC vector; plot_fsc')
        fname_eps  = trim(tmpl_fname)//'.eps'
        fname_pdf  = trim(tmpl_fname)//'.pdf'
        call CPlot2D__new(plot2D, trim(tmpl_fname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.d0)
        call CPlot2D__SetYAxisSize(plot2D, 400.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_FALSE)
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,0.d0,1.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet, 1.0/res(k), fsc(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        title = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%to_char())
        title = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%to_char())
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(ps2pdf_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_fsc

    subroutine plot_fsc2( n, fsc1, fsc2, res, smpd, tmpl_fname )
        integer,           intent(in) :: n
        real,              intent(in) :: fsc1(n), fsc2(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        type(string)              :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet1
        type(CDataSet_type)       :: dataSet2
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, fname_eps
        integer  :: k,iostat
        if( n == 0 ) THROW_HARD('Empty FSC vector; plot_fsc')
        fname_eps  = trim(tmpl_fname)//'.eps'
        fname_pdf  = trim(tmpl_fname)//'.pdf'
        call CPlot2D__new(plot2D, trim(tmpl_fname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.d0)
        call CPlot2D__SetYAxisSize(plot2D, 400.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_FALSE)
        call CDataSet__new(dataSet1)
        do k = 1,n
            call CDataSet_addpoint(dataSet1, 1.0/res(k), fsc1(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataSet1)
        call CDataSet__delete(dataSet1)
        call CDataSet__new(dataSet2)
        call CDataSet__SetDatasetColor(dataSet2, 0.d0,0.d0,1.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet2, 1.0/res(k), fsc2(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataSet2)
        call CDataSet__delete(dataSet2)
        title = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%to_char())
        title = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%to_char())
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(ps2pdf_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_fsc2

    subroutine plot_phrand_fsc( n, fsc, fsc_t, fsc_n, res, smpd, tmpl_fname )
        integer,           intent(in) :: n
        real,              intent(in) :: fsc(n), fsc_n(n), fsc_t(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        type(string)              :: title
        type(CPlot2D_type)        :: plot2D
        type(CDataSet_type)       :: dataSet
        character(len=LONGSTRLEN) :: ps2pdf_cmd, fname_pdf, fname_eps
        integer  :: k,iostat
        if( n == 0 ) THROW_HARD('Empty FSC vector; plot_fsc')
        fname_eps  = trim(tmpl_fname)//'.eps'
        fname_pdf  = trim(tmpl_fname)//'.pdf'
        call CPlot2D__new(plot2D, trim(tmpl_fname)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 400.d0)
        call CPlot2D__SetYAxisSize(plot2D, 400.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetFlipY(plot2D, C_FALSE)
        ! corrected fsc
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,0.d0,1.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet, 1.0/res(k), fsc(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        ! raw fsc
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 0.d0,1.d0,0.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet, 1.0/res(k), fsc_t(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        ! phase randomized fsc
        call CDataSet__new(dataSet)
        call CDataSet__SetDrawMarker(dataSet, C_FALSE)
        call CDataSet__SetDatasetColor(dataSet, 1.d0,0.d0,0.d0)
        do k = 1,n
            call CDataSet_addpoint(dataSet, 1.0/res(k), fsc_n(k))
        end do
        call CPlot2D__AddDataSet(plot2D, dataset)
        call CDataSet__delete(dataset)
        title = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%to_char())
        title = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%to_char())
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(ps2pdf_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_phrand_fsc

    subroutine plot_fsc_area_score( result, fbody )
        type(fsc_area_score_result), intent(in) :: result
        character(len=*),            intent(in) :: fbody
        type(string)              :: title
        type(CPlot2D_type)        :: plot2D
        character(len=LONGSTRLEN) :: fname_eps, fname_png, plot_title
        character(len=XLONGSTRLEN):: ps2png_cmd
        real, allocatable :: freq(:), mean(:), stdev(:), minv(:), maxv(:), lower(:), upper(:), crossings(:)
        logical, allocatable :: valid_shell(:)
        real(dp) :: sumv, sumsq, mu, sig
        real     :: val, max_freq, best_res, worst_res
        integer  :: n, k, idir, cnt, ibin, ncross, iostat
        if( result%nfreqs <= 0 ) THROW_HARD('Empty FSC area score result; plot_fsc_area_score')
        if( .not.allocated(result%res) .or. .not.allocated(result%cfsc) .or. .not.allocated(result%counts) )then
            THROW_HARD('Incomplete FSC area score result; plot_fsc_area_score')
        endif
        n = result%nfreqs
        allocate(freq(n), mean(n), stdev(n), minv(n), maxv(n), lower(n), upper(n), crossings(n), source=0.)
        allocate(valid_shell(n), source=.false.)
        do k = 1,n
            if( result%res(k) > 0. ) freq(k) = 1. / result%res(k)
            sumv = 0.d0
            sumsq = 0.d0
            cnt = 0
            minv(k) = huge(1.)
            maxv(k) = -huge(1.)
            do idir = 1,result%ndirs
                if( result%counts(k,idir) < result%min_count ) cycle
                val = max(0., min(1., result%cfsc(k,idir)))
                sumv = sumv + real(val,dp)
                sumsq = sumsq + real(val,dp) * real(val,dp)
                cnt = cnt + 1
                minv(k) = min(minv(k), val)
                maxv(k) = max(maxv(k), val)
            end do
            if( cnt > 0 )then
                mu = sumv / real(cnt,dp)
                sig = sqrt(max(0.d0, sumsq / real(cnt,dp) - mu * mu))
                mean(k)  = real(mu)
                stdev(k) = real(sig)
                valid_shell(k) = .true.
            else
                mean(k)  = 0.
                stdev(k) = 0.
                minv(k)  = 0.
                maxv(k)  = 0.
            endif
            lower(k) = max(0., mean(k) - stdev(k))
            upper(k) = min(1., mean(k) + stdev(k))
        end do
        ncross   = 0
        best_res = huge(1.)
        worst_res = -huge(1.)
        if( allocated(result%crossing_find) .and. allocated(result%crossing_res) )then
            do idir = 1,result%ndirs
                if( result%crossing_find(idir) <= 0. .or. result%crossing_res(idir) <= 0. ) cycle
                ibin = nint(result%crossing_find(idir))
                if( ibin < 1 .or. ibin > n ) cycle
                crossings(ibin) = crossings(ibin) + 1.
                ncross = ncross + 1
                best_res  = min(best_res,  result%crossing_res(idir))
                worst_res = max(worst_res, result%crossing_res(idir))
            end do
            if( ncross > 0 ) crossings = crossings / real(ncross)
        endif
        max_freq = maxval(freq)
        if( max_freq <= 0. ) THROW_HARD('Invalid resolution axis; plot_fsc_area_score')
        if( ncross > 0 )then
            write(plot_title,'(A,F6.3,A,F6.2,A,F6.2,A)') 'Directional FSC area score  cFAR=', &
                &result%cfar, '  worst=', worst_res, ' A  best=', best_res, ' A'
        else
            write(plot_title,'(A,F6.3)') 'Directional FSC area score  cFAR=', result%cfar
        endif
        fname_eps = trim(fbody)//'_plot.eps'
        fname_png = trim(fbody)//'_plot.png'
        call CPlot2D__new(plot2D, trim(plot_title)//C_NULL_CHAR)
        call CPlot2D__SetXAxisSize(plot2D, 560.d0)
        call CPlot2D__SetYAxisSize(plot2D, 320.d0)
        call CPlot2D__SetDrawLegend(plot2D, C_FALSE)
        call CPlot2D__SetDrawXAxisGridLines(plot2D, C_FALSE)
        call CPlot2D__SetDrawYAxisGridLines(plot2D, C_TRUE)
        call CPlot2D__SetFlipY(plot2D, C_FALSE)
        call add_invisible_bounds
        call add_filled_band(minv, maxv, 0.71_c_double, 0.83_c_double, 0.96_c_double)
        call add_filled_band(lower, upper, 0.54_c_double, 0.73_c_double, 0.91_c_double)
        call add_crossing_bars
        call add_segment(0., result%threshold, max_freq, result%threshold, &
            &0.53_c_double, 0.53_c_double, 0.50_c_double, 1.0_c_double, C_TRUE)
        call add_curve(mean, 0.09_c_double, 0.37_c_double, 0.65_c_double, 2.2_c_double, C_FALSE)
        title = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%to_char())
        title = 'FSC / relative occurrence'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%to_char())
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ps2png_cmd = 'gs -q -sDEVICE=png16m -dNOPAUSE -dBATCH -dSAFER -dEPSCrop -r144 '&
            &//'-dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sOutputFile='//trim(fname_png)//' '//trim(fname_eps)
        call exec_cmdline(ps2png_cmd, suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
        write(logfhandle,'(A,1X,A)') '>>> FSC area score plot:', trim(fname_png)
        deallocate(freq, mean, stdev, minv, maxv, lower, upper, crossings, valid_shell)
    contains
        subroutine add_invisible_bounds
            type(CDataSet_type) :: dataSet
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDrawLine(dataSet, C_FALSE)
            call CDataSet_addpoint(dataSet, 0., 0.)
            call CDataSet_addpoint(dataSet, max_freq, 1.)
            call CPlot2D__AddDataSet(plot2D, dataSet)
            call CDataSet__delete(dataSet)
        end subroutine add_invisible_bounds

        subroutine add_filled_band( lo, hi, r, g, b )
            real,           intent(in) :: lo(n), hi(n)
            real(C_double), intent(in) :: r, g, b
            integer :: first, last, kk
            first = 0
            last  = 0
            do kk = 1,n
                if( valid_shell(kk) )then
                    if( first == 0 ) first = kk
                    last = kk
                else if( first > 0 )then
                    call add_filled_band_segment(lo, hi, first, last, r, g, b)
                    first = 0
                    last  = 0
                endif
            end do
            if( first > 0 ) call add_filled_band_segment(lo, hi, first, last, r, g, b)
        end subroutine add_filled_band

        subroutine add_filled_band_segment( lo, hi, first, last, r, g, b )
            real,           intent(in) :: lo(n), hi(n)
            integer,        intent(in) :: first, last
            real(C_double), intent(in) :: r, g, b
            type(CDataSet_type) :: dataSet
            integer :: kk
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDrawLine(dataSet, C_FALSE)
            call CDataSet__SetFillArea(dataSet, C_TRUE)
            call CDataSet__SetDatasetColor(dataSet, r, g, b)
            do kk = first,last
                call CDataSet_addpoint(dataSet, freq(kk), hi(kk))
            end do
            do kk = last,first,-1
                call CDataSet_addpoint(dataSet, freq(kk), lo(kk))
            end do
            call CPlot2D__AddDataSet(plot2D, dataSet)
            call CDataSet__delete(dataSet)
        end subroutine add_filled_band_segment

        subroutine add_crossing_bars
            real :: left_edge, right_edge, half_width
            integer :: kk
            if( ncross == 0 ) return
            do kk = 1,n
                if( crossings(kk) <= 0. ) cycle
                call bar_edges(kk, left_edge, right_edge)
                half_width = 0.43 * (right_edge - left_edge)
                call add_filled_rect(max(0., freq(kk) - half_width), freq(kk) + half_width, 0., crossings(kk), &
                    &0.23_c_double, 0.43_c_double, 0.07_c_double)
            end do
        end subroutine add_crossing_bars

        subroutine bar_edges( kk, left_edge, right_edge )
            integer, intent(in)  :: kk
            real,    intent(out) :: left_edge, right_edge
            if( n == 1 )then
                left_edge  = max(0., freq(kk) - 0.01)
                right_edge = freq(kk) + 0.01
            else if( kk == 1 )then
                left_edge  = max(0., freq(kk) - 0.5 * (freq(kk+1) - freq(kk)))
                right_edge = 0.5 * (freq(kk) + freq(kk+1))
            else if( kk == n )then
                left_edge  = 0.5 * (freq(kk-1) + freq(kk))
                right_edge = freq(kk) + 0.5 * (freq(kk) - freq(kk-1))
            else
                left_edge  = 0.5 * (freq(kk-1) + freq(kk))
                right_edge = 0.5 * (freq(kk) + freq(kk+1))
            endif
        end subroutine bar_edges

        subroutine add_filled_rect( xlo, xhi, ylo, yhi, r, g, b )
            real,           intent(in) :: xlo, xhi, ylo, yhi
            real(C_double), intent(in) :: r, g, b
            type(CDataSet_type) :: dataSet
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDrawLine(dataSet, C_FALSE)
            call CDataSet__SetFillArea(dataSet, C_TRUE)
            call CDataSet__SetDatasetColor(dataSet, r, g, b)
            call CDataSet_addpoint(dataSet, xlo, ylo)
            call CDataSet_addpoint(dataSet, xlo, yhi)
            call CDataSet_addpoint(dataSet, xhi, yhi)
            call CDataSet_addpoint(dataSet, xhi, ylo)
            call CPlot2D__AddDataSet(plot2D, dataSet)
            call CDataSet__delete(dataSet)
        end subroutine add_filled_rect

        subroutine add_curve( y, r, g, b, linewidth, dashed )
            real,           intent(in) :: y(n)
            real(C_double), intent(in) :: r, g, b, linewidth
            logical(C_bool), intent(in):: dashed
            integer :: first, last, kk
            first = 0
            last  = 0
            do kk = 1,n
                if( valid_shell(kk) )then
                    if( first == 0 ) first = kk
                    last = kk
                else if( first > 0 )then
                    call add_curve_segment(y, first, last, r, g, b, linewidth, dashed)
                    first = 0
                    last  = 0
                endif
            end do
            if( first > 0 ) call add_curve_segment(y, first, last, r, g, b, linewidth, dashed)
        end subroutine add_curve

        subroutine add_curve_segment( y, first, last, r, g, b, linewidth, dashed )
            real,           intent(in) :: y(n)
            integer,        intent(in) :: first, last
            real(C_double), intent(in) :: r, g, b, linewidth
            logical(C_bool), intent(in):: dashed
            type(CDataSet_type) :: dataSet
            integer :: kk
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDrawLine(dataSet, C_TRUE)
            call CDataSet__SetLineWidth(dataSet, linewidth)
            call CDataSet__SetDashedLine(dataSet, dashed)
            call CDataSet__SetDatasetColor(dataSet, r, g, b)
            do kk = first,last
                call CDataSet_addpoint(dataSet, freq(kk), y(kk))
            end do
            call CPlot2D__AddDataSet(plot2D, dataSet)
            call CDataSet__delete(dataSet)
        end subroutine add_curve_segment

        subroutine add_segment( x1, y1, x2, y2, r, g, b, linewidth, dashed )
            real,           intent(in) :: x1, y1, x2, y2
            real(C_double), intent(in) :: r, g, b, linewidth
            logical(C_bool), intent(in):: dashed
            type(CDataSet_type) :: dataSet
            call CDataSet__new(dataSet)
            call CDataSet__SetDrawMarker(dataSet, C_FALSE)
            call CDataSet__SetDrawLine(dataSet, C_TRUE)
            call CDataSet__SetLineWidth(dataSet, linewidth)
            call CDataSet__SetDashedLine(dataSet, dashed)
            call CDataSet__SetDatasetColor(dataSet, r, g, b)
            call CDataSet_addpoint(dataSet, x1, y1)
            call CDataSet_addpoint(dataSet, x2, y2)
            call CPlot2D__AddDataSet(plot2D, dataSet)
            call CDataSet__delete(dataSet)
        end subroutine add_segment
    end subroutine plot_fsc_area_score

    subroutine fibonacci_sphere_dirs( ndirs, dirs )
        integer, intent(in)  :: ndirs
        real,    intent(out) :: dirs(3,ndirs)
        real(dp), parameter :: GOLDEN_ANGLE = DPI * (3.d0 - sqrt(5.d0))
        real(dp) :: z, rho, phi
        integer  :: idir
        do idir = 1, ndirs
            z          = 1.d0 - 2.d0 * (real(idir,dp) - 0.5d0) / real(ndirs,dp)
            rho        = sqrt(max(0.d0, 1.d0 - z*z))
            phi        = real(idir - 1, dp) * GOLDEN_ANGLE
            dirs(1,idir) = real(rho * cos(phi))
            dirs(2,idir) = real(rho * sin(phi))
            dirs(3,idir) = real(z)
        end do
    end subroutine fibonacci_sphere_dirs

    subroutine score_cfsc_curve( n, cfsc, counts, res, threshold, min_count, wauc, crossing_find, crossing_res, included_bins )
        integer, intent(in)  :: n, min_count
        real,    intent(in)  :: cfsc(n), res(n), threshold
        integer, intent(in)  :: counts(n)
        real,    intent(out) :: wauc, crossing_find, crossing_res
        integer, intent(out) :: included_bins
        real(dp) :: area
        integer  :: k
        area          = 0.d0
        crossing_find = 0.
        crossing_res  = 0.
        included_bins = 0
        do k = 1, n
            if( counts(k) < min_count ) cycle
            if( cfsc(k) <= threshold )then
                crossing_find = real(k)
                crossing_res  = res(k)
                exit
            endif
            area = area + real(cfsc(k),dp) * real(k*k,dp)
            included_bins = included_bins + 1
        end do
        wauc = real(area)
    end subroutine score_cfsc_curve

    subroutine write_fsc_area_score_outputs( result, fbody )
        type(fsc_area_score_result), intent(in) :: result
        character(len=*),            intent(in) :: fbody
        character(len=LONGSTRLEN) :: summary_fname, curves_fname, dirs_fname
        integer :: funit, iostat, idir, k
        summary_fname = trim(fbody)//'_summary.txt'
        curves_fname  = trim(fbody)//'_curves.csv'
        dirs_fname    = trim(fbody)//'_directions.csv'
        call fopen(funit, file=string(summary_fname), status='REPLACE', action='WRITE', iostat=iostat)
        call fileiochk('write_fsc_area_score_outputs; open '//trim(summary_fname), iostat)
        write(funit,'(A)') 'FSC AREA SCORE / CONICAL FSC AREA RATIO'
        write(funit,'(A,1X,I0)')      'direction_axes', result%ndirs
        write(funit,'(A,1X,I0)')      'frequency_shells', result%nfreqs
        write(funit,'(A,1X,F10.4)')   'cone_half_angle_deg', result%cone_half_angle_deg
        write(funit,'(A,1X,F10.4)')   'threshold', result%threshold
        write(funit,'(A,1X,I0)')      'min_voxels_per_shell_cone', result%min_count
        write(funit,'(A,1X,F12.6)')   'cfar', result%cfar
        write(funit,'(A,1X,ES16.8)')  'wauc_min', result%wauc_min
        write(funit,'(A,1X,ES16.8)')  'wauc_max', result%wauc_max
        write(funit,'(A,1X,I0)')      'min_direction', result%min_dir
        write(funit,'(A,1X,I0)')      'max_direction', result%max_dir
        call fclose(funit)
        call fopen(funit, file=string(curves_fname), status='REPLACE', action='WRITE', iostat=iostat)
        call fileiochk('write_fsc_area_score_outputs; open '//trim(curves_fname), iostat)
        write(funit,'(A)', advance='no') 'wave_number,resolution_A'
        do idir = 1, result%ndirs
            write(funit,'(A)', advance='no') ',dir'//int2str_pad(idir,3)
        end do
        write(funit,'(A)') ''
        do k = 1, result%nfreqs
            write(funit,'(I0,A,F12.5)', advance='no') k, ',', result%res(k)
            do idir = 1, result%ndirs
                if( result%counts(k,idir) >= result%min_count )then
                    write(funit,'(A,ES16.8)', advance='no') ',', result%cfsc(k,idir)
                else
                    write(funit,'(A)', advance='no') ',nan'
                endif
            end do
            write(funit,'(A)') ''
        end do
        call fclose(funit)
        call fopen(funit, file=string(dirs_fname), status='REPLACE', action='WRITE', iostat=iostat)
        call fileiochk('write_fsc_area_score_outputs; open '//trim(dirs_fname), iostat)
        write(funit,'(A)') 'direction,x,y,z,wauc,crossing_wave_number,crossing_resolution_A,included_bins'
        do idir = 1, result%ndirs
            write(funit,'(I0,A,F12.7,A,F12.7,A,F12.7,A,ES16.8,A,F10.3,A,F12.5,A,I0)') &
                &idir, ',', result%dirs(1,idir), ',', result%dirs(2,idir), ',', result%dirs(3,idir), ',', &
                &result%wauc(idir), ',', result%crossing_find(idir), ',', result%crossing_res(idir), ',', &
                &result%included_bins(idir)
        end do
        call fclose(funit)
        write(logfhandle,'(A,1X,A)') '>>> FSC area score summary:', trim(summary_fname)
        write(logfhandle,'(A,1X,A)') '>>> FSC area score curves:',  trim(curves_fname)
        write(logfhandle,'(A,1X,A)') '>>> FSC area score directions:', trim(dirs_fname)
    end subroutine write_fsc_area_score_outputs

    ! TYPE fsc_area_score_result BOUND PROCEDURES

    subroutine new( self, vol, ndirs, cone_half_angle_deg, threshold, min_count)
        class(fsc_area_score_result), intent(out) :: self
        class(image),                 intent(in)  :: vol
        integer,                      intent(in)  :: ndirs
        real,                         intent(in)  :: cone_half_angle_deg, threshold
        integer,                      intent(in)  :: min_count
        if( ndirs < 1 ) THROW_HARD('number of directions must be positive; calc_fsc_area_score')
        if( cone_half_angle_deg <= 0. .or. cone_half_angle_deg >= 90. )then
            THROW_HARD('cone half-angle must be in (0,90) degrees; calc_fsc_area_score')
        endif
        if( threshold < -1. .or. threshold > 1. ) THROW_HARD('threshold outside [-1,1]; calc_fsc_area_score')
        call self%kill
        self%nfreqs              = vol%get_filtsz()
        self%ndirs               = ndirs
        self%min_count           = max(1, min_count)
        self%cone_half_angle_deg = cone_half_angle_deg
        self%threshold           = threshold
        allocate(self%res(self%nfreqs), self%dirs(3,self%ndirs), self%cfsc(self%nfreqs,self%ndirs), &
        &self%counts(self%nfreqs,self%ndirs), self%wauc(self%ndirs), self%crossing_find(self%ndirs),&
        &self%crossing_res(self%ndirs), self%included_bins(self%ndirs))
        self%res = vol%get_res()
        call fibonacci_sphere_dirs(self%ndirs, self%dirs)
        self%exists = .true.
    end subroutine new

    subroutine calc_fsc_area_score( self, even, odd, state )
        class(fsc_area_score_result), intent(inout) :: self
        class(image),                 intent(inout) :: even, odd
        integer,            optional, intent(in)    :: state
        logical, parameter :: VERBOSE = .false.
        integer :: n, idir, loc(1)
        if( .not. even%is_ft() ) call even%fft()
        if( .not. odd%is_ft()  ) call odd%fft()
        call even%conical_fsc(odd, self%dirs, self%cone_half_angle_deg, self%min_count, self%cfsc, self%counts)
        do idir = 1, self%ndirs
            call score_cfsc_curve(self%nfreqs, self%cfsc(:,idir), self%counts(:,idir), self%res, self%threshold, self%min_count, &
                &self%wauc(idir), self%crossing_find(idir), self%crossing_res(idir), self%included_bins(idir))
        end do
        loc = minloc(self%wauc)
        self%min_dir  = loc(1)
        self%wauc_min = self%wauc(self%min_dir)
        loc = maxloc(self%wauc)
        self%max_dir  = loc(1)
        self%wauc_max = self%wauc(self%max_dir)
        if( self%wauc_max > 0. )then
            self%cfar = self%wauc_min / self%wauc_max
        endif
        if( present(state) )then
            write(logfhandle,'(A,I2,A,1X,F8.4)') '>>> FSC AREA SCORE / CONICAL FSC AREA RATIO (cFAR) STATE ', state,':',self%cfar
        else
            write(logfhandle,'(A,1X,F8.4)') '>>> FSC AREA SCORE / CONICAL FSC AREA RATIO (cFAR):', self%cfar
        endif
        if( VERBOSE )then
            write(logfhandle,'(A,1X,I0)')      '    Direction axes:', self%ndirs
            write(logfhandle,'(A,1X,F7.3)')    '    Cone half-angle (degrees):', self%cone_half_angle_deg
            write(logfhandle,'(A,1X,F7.3)')    '    FSC threshold:', self%threshold
            write(logfhandle,'(A,1X,ES12.5)')  '    Minimum weighted area:', self%wauc_min
            write(logfhandle,'(A,1X,ES12.5)')  '    Maximum weighted area:', self%wauc_max
        endif
    end subroutine calc_fsc_area_score

    subroutine kill( self )
        class(fsc_area_score_result), intent(out) :: self
        self%nfreqs              = 0
        self%ndirs               = 0
        self%min_count           = 1
        self%min_dir             = 0
        self%max_dir             = 0
        self%cone_half_angle_deg = 20.
        self%threshold           = 0.143
        self%cfar                = 0.
        self%wauc_min            = 0.
        self%wauc_max            = 0.
        if( allocated(self%res) )then
            deallocate(self%res, self%dirs, self%cfsc, self%wauc, self%crossing_find, &
            &self%crossing_res, self%counts, self%included_bins)
        endif
        self%exists = .false.
    end subroutine kill

end module simple_fsc
