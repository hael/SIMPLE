!@descr: various Fourier Shell Correlation utilities
module simple_fsc
use simple_core_module_api
use simple_image, only: image
use CPlot2D_wrapper_module
implicit none

public :: phase_rand_fsc, plot_fsc, plot_fsc2, plot_phrand_fsc
public :: fsc_area_score_result, calc_fsc_area_score, write_fsc_area_score_outputs
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

    subroutine calc_fsc_area_score( even, odd, ndirs, cone_half_angle_deg, threshold, min_count, result )
        class(image),                intent(inout) :: even, odd
        integer,                     intent(in)    :: ndirs
        real,                        intent(in)    :: cone_half_angle_deg, threshold
        integer,                     intent(in)    :: min_count
        type(fsc_area_score_result), intent(out)   :: result
        integer :: n, idir, loc(1), min_count_eff
        if( ndirs < 1 ) THROW_HARD('number of directions must be positive; calc_fsc_area_score')
        if( cone_half_angle_deg <= 0. .or. cone_half_angle_deg >= 90. )then
            THROW_HARD('cone half-angle must be in (0,90) degrees; calc_fsc_area_score')
        endif
        if( threshold < -1. .or. threshold > 1. ) THROW_HARD('threshold outside [-1,1]; calc_fsc_area_score')
        if( .not. even%is_ft() ) call even%fft()
        if( .not. odd%is_ft()  ) call odd%fft()
        n             = even%get_filtsz()
        min_count_eff = max(1, min_count)
        result%nfreqs              = n
        result%ndirs               = ndirs
        result%min_count           = min_count_eff
        result%cone_half_angle_deg = cone_half_angle_deg
        result%threshold           = threshold
        allocate(result%res(n), result%dirs(3,ndirs), result%cfsc(n,ndirs), result%counts(n,ndirs))
        allocate(result%wauc(ndirs), result%crossing_find(ndirs), result%crossing_res(ndirs), result%included_bins(ndirs))
        result%res = even%get_res()
        call fibonacci_sphere_dirs(ndirs, result%dirs)
        call even%conical_fsc(odd, result%dirs, cone_half_angle_deg, min_count_eff, result%cfsc, result%counts)
        do idir = 1, ndirs
            call score_cfsc_curve(n, result%cfsc(:,idir), result%counts(:,idir), result%res, threshold, min_count_eff, &
                &result%wauc(idir), result%crossing_find(idir), result%crossing_res(idir), result%included_bins(idir))
        end do
        loc = minloc(result%wauc)
        result%min_dir  = loc(1)
        result%wauc_min = result%wauc(result%min_dir)
        loc = maxloc(result%wauc)
        result%max_dir  = loc(1)
        result%wauc_max = result%wauc(result%max_dir)
        if( result%wauc_max > 0. )then
            result%cfar = result%wauc_min / result%wauc_max
        else
            result%cfar = 0.
        endif
        write(logfhandle,'(A)') '>>> FSC AREA SCORE / CONICAL FSC AREA RATIO'
        write(logfhandle,'(A,1X,I0)')      '    Direction axes:', ndirs
        write(logfhandle,'(A,1X,F7.3)')    '    Cone half-angle (degrees):', cone_half_angle_deg
        write(logfhandle,'(A,1X,F7.3)')    '    FSC threshold:', threshold
        write(logfhandle,'(A,1X,ES12.5)')  '    Minimum weighted area:', result%wauc_min
        write(logfhandle,'(A,1X,ES12.5)')  '    Maximum weighted area:', result%wauc_max
        write(logfhandle,'(A,1X,F8.4)')    '    cFAR:', result%cfar
    end subroutine calc_fsc_area_score

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

end module simple_fsc
