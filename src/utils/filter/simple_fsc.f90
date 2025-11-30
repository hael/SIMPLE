module simple_fsc
include 'simple_lib.f08'
use simple_image, only: image
use CPlot2D_wrapper_module
implicit none

public :: phase_rand_fsc, plot_fsc, plot_fsc2, plot_phrand_fsc
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

end module simple_fsc
