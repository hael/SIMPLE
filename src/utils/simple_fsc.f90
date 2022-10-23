module simple_fsc
include 'simple_lib.f08'
use simple_image, only: image
use CPlot2D_wrapper_module
implicit none

public :: phase_rand_fsc, plot_fsc, plot_phrand_fsc, phaseplate_correct_fsc
private
#include "simple_local_flags.inc"

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
        call even%zero_background
        call odd%zero_background
        call even%mul(envmask)                  ! mask
        call odd%mul(envmask)
        call even%fft()                         ! Fourier space
        call odd%fft()
        call even%fsc(odd, fsc_t)               ! FSC
        ! Randomize then calculate masked FSC
        k_rand = get_lplim_at_corr(fsc_t, ENVMSK_FSC_THRESH)
        if( k_rand > n-3 )then
            ! reverts to sherical masking
            call even%ifft()
            call odd%ifft()
            call even%mask(msk, 'soft')
            call odd%mask(msk,  'soft')
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
            call even%zero_background
            call odd%zero_background
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
            call arr2file(fsc_t, 'fsct_state'//int2str_pad(state,2)//BIN_EXT)
            call arr2file(fsc_n, 'fscn_state'//int2str_pad(state,2)//BIN_EXT)
        endif
    end subroutine phase_rand_fsc

    subroutine plot_fsc( n, fsc, res, smpd, tmpl_fname )
        integer,           intent(in) :: n
        real,              intent(in) :: fsc(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        real, parameter           :: SCALE = 40.
        type(str4arr)             :: title
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
        title%str = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_fsc

    subroutine plot_phrand_fsc( n, fsc, fsc_t, fsc_n, res, smpd, tmpl_fname )
        integer,           intent(in) :: n
        real,              intent(in) :: fsc(n), fsc_n(n), fsc_t(n), res(n), smpd
        character(len=*),  intent(in) :: tmpl_fname
        real, parameter           :: SCALE = 40.
        type(str4arr)             :: title
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
        title%str = 'Resolution (Angstroms^-1)'//C_NULL_CHAR
        call CPlot2D__SetXAxisTitle(plot2D, title%str)
        title%str = 'Fourier Shell Correlations'//C_NULL_CHAR
        call CPlot2D__SetYAxisTitle(plot2D, title%str)
        call CPlot2D__OutputPostScriptPlot(plot2D, trim(fname_eps)//C_NULL_CHAR)
        call CPlot2D__delete(plot2D)
        ! conversion to PDF
        ps2pdf_cmd = 'gs -q -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600 -sOutputFile='&
            &//trim(fname_pdf)//' '//trim(fname_eps)
        call exec_cmdline(trim(adjustl(ps2pdf_cmd)), suppress_errors=.true., exitstat=iostat)
        if( iostat == 0 ) call del_file(fname_eps)
    end subroutine plot_phrand_fsc

    subroutine phaseplate_correct_fsc( fsc, find_plate )
        real,    intent(inout) :: fsc(:)
        integer, intent(out)   :: find_plate
        logical, allocatable :: peakpos(:)
        integer :: k, n
        real    :: peakavg
        n = size(fsc)
        ! find FSC peaks
        peakpos = peakfinder(fsc)
        ! filter out all peaks FSC < 0.5
        where(fsc < 0.5) peakpos = .false.
        ! calculate peak average
        peakavg = sum(fsc, mask=peakpos)/real(count(peakpos))
        ! identify peak with highest frequency
        do k=n,1,-1
            if( peakpos(k) )then
                find_plate = k
                exit
            endif
        end do
        ! replace with average FSC peak value up to last peak (find_plate)
        do k=1,find_plate
            fsc(k) = peakavg
        end do
    end subroutine phaseplate_correct_fsc

end module simple_fsc
