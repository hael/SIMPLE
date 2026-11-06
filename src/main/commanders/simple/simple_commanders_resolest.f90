!@descr: for resolution estimation
module simple_commanders_resolest
use simple_commanders_api
use simple_pftc_srch_api
use simple_fsc
use simple_commanders_euclid, only: commander_calc_pspec
use simple_refine3D_fnames,  only: refine3D_fsc_fbody, refine3D_fsc_fname
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_fsc
  contains
    procedure :: execute      => exec_fsc
end type commander_fsc

type, extends(commander_base) :: commander_fsc_area_score
  contains
    procedure :: execute      => exec_fsc_area_score
end type commander_fsc_area_score

type, extends(commander_base) :: commander_uniform_filter2D
  contains
    procedure :: execute      => exec_uniform_filter2D
end type commander_uniform_filter2D

type, extends(commander_base) :: commander_uniform_filter3D
  contains
    procedure :: execute      => exec_uniform_filter3D
end type commander_uniform_filter3D

type, extends(commander_base) :: commander_nu_filt3D
    contains
        procedure :: execute      => exec_nu_filt3D
end type commander_nu_filt3D

type, extends(commander_base) :: commander_dose_filter3D
  contains
    procedure :: execute      => exec_dose_filter3D
end type commander_dose_filter3D

type, extends(commander_base) :: commander_estimate_lpstages
  contains
    procedure :: execute      => exec_estimate_lpstages
end type commander_estimate_lpstages

contains

    !> accumulates a dose series of half-map pairs under a measured attenuation model
    subroutine exec_dose_filter3D( self, cline )
        use simple_dose_series,   only: dose_series, shell_independence, shell_correlation, DOSE_WEIGHTS_MODEL, DOSE_WEIGHTS_CONTROL, &
                                      &DOSE_WEIGHTS_FLAT
        use simple_dose_series_manifest, only: dose_series_manifest
        use simple_dose_curve,    only: dose_curve
        use simple_estimate_ssnr, only: get_resolution
        use simple_syslib,        only: simple_abspath
        use simple_string_utils,  only: posix_quote, fortran_quote
        class(commander_dose_filter3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        !> one window of a dose series, as the manifest describes it
        type :: dose_window
            type(string) :: even, odd
            integer      :: ffrom = 0, fto = 0
            real(dp)     :: dose_lo = 0._dp, dose_hi = 0._dp
        end type dose_window
        type(parameters)   :: params
        type(dose_series)  :: series
        type(dose_curve)   :: curve
        type(dose_series_manifest) :: input_series
        type(dose_window), allocatable :: win(:)
        complex(sp), allocatable :: em(:,:), om(:,:), ae(:), ao(:), s0(:)
        integer,     allocatable :: shell(:)
        real(dp),    allocatable :: kfreq(:)
        real,        allocatable :: fsc(:), fsc_msk(:), res(:), fsc_tab(:,:)
        type(image)  :: vol
        type(image)  :: mskimg
        type(string) :: outbody, manifest, message
        real(dp)     :: mask_frac, shell_corr, hdr_smpd, hdr_kv
        real         :: fsc05(4), fsc0143(4), mfsc05(4), mfsc0143(4)
        integer      :: hdr_box, nwin, w, nsh, ncoeff, lp_shell, hp_shell, klo, khi
        ! the resolutions of the weights plot, in Angstroms
        real(dp), parameter :: PLOT_RES(5) = [20._dp, 10._dp, 8._dp, 6._dp, 4._dp]
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('bfac_model') ) call cline%set('bfac_model', 'line')
        if( .not. cline%defined('outvol')     ) call cline%set('outvol', 'dose_series.mrc')
        ! B-factors carry no information below ~20 A, where the mask and solvent artefacts live
        if( .not. cline%defined('hp')         ) call cline%set('hp', 20.)
        if( .not. cline%defined('filetab')    ) THROW_HARD('dose_filter3D needs filetab, the dose series manifest')
        ! the manifest, not a volume key, carries the box: parse it before the parameters are made;
        ! absolute, because params%new enters a fresh run directory
        manifest = cline%get_carg('filetab')
        if( .not. file_exists(manifest) ) THROW_HARD('Dose series manifest not found')
        manifest = simple_abspath(manifest)
        if( .not. input_series%read(manifest, message) ) THROW_HARD(message%to_char())
        if( .not. input_series%select_windows(message) ) THROW_HARD(message%to_char())
        call input_series%get_geometry(hdr_box, hdr_smpd, hdr_kv)
        allocate(win(input_series%get_nwindows()))
        do w = 1, size(win)
            call input_series%get_window(w, win(w)%ffrom, win(w)%fto, win(w)%dose_lo, win(w)%dose_hi, &
                &win(w)%even, win(w)%odd)
            write(logfhandle,'(A,2(I0,1X),A,2(ES14.6,1X))') '>>> DOSE SERIES: FRACTIONS ', &
                &win(w)%ffrom, win(w)%fto, 'PHYSICAL DOSE BOUNDS ', win(w)%dose_lo, win(w)%dose_hi
        enddo
        call input_series%kill
        call cline%set('box', hdr_box)
        call params%new(cline)
        if( trim(params%bfac_model).ne.'line' .and. trim(params%bfac_model).ne.'free' )then
            THROW_HARD('bfac_model must be line or free')
        endif
        if( abs(hdr_smpd - real(params%smpd,dp)) > 1.e-4_dp ) THROW_HARD('Manifest sampling disagrees with smpd')
        nwin = size(win)
        write(logfhandle,'(A,I0,A)') '>>> DOSE SERIES: ', nwin, ' WINDOWS'
        ! the mask the cross-spectra and the cross-validation see; its soft edge extends past msk
        if( params%msk + COSMSKHALFWIDTH > real(params%box/2) ) THROW_HARD('mskdiam does not fit the box')
        ! the indicator of that mask, for the independent count; backgr=0 so nothing is subtracted first
        call mskimg%new([params%box,params%box,params%box], params%smpd)
        mskimg = 1.
        call mskimg%mask3D_soft(params%msk, backgr=0.)
        mask_frac  = shell_independence(mskimg%get_rmat())
        if( mask_frac <= 0._dp ) THROW_HARD('Mask leaves no support')
        shell_corr = shell_correlation(mskimg%get_rmat())
        write(logfhandle,'(A,F7.4)') '>>> DOSE SERIES: INDEPENDENT FRACTION OF A SHELL UNDER THE MASK ', mask_frac
        write(logfhandle,'(A,F7.4)') '>>> DOSE SERIES: VARIANCE INFLATION ACROSS SHELLS UNDER THE MASK ', shell_corr
        call vol%new([params%box,params%box,params%box], params%smpd)
        nsh = vol%get_filtsz()
        call lattice_shells(vol, nsh, shell, ncoeff)
        allocate(kfreq(nsh))
        do w = 1, nsh
            kfreq(w) = real(w,dp) / (real(params%box,dp) * real(params%smpd,dp))
        enddo
        call series%new(nwin, kfreq, hdr_kv, trim(params%bfac_model).eq.'line')
        do w = 1, nwin
            call series%set_window(w, win(w)%dose_lo, win(w)%dose_hi)
        enddo
        ! pass one: the masked halves, for the statistics; they stay until the paired cross-validation
        call gather_series(win, params, ncoeff, .true., em, om)
        call series%accumulate_gram(em, om, shell, mask_frac, shell_corr)
        lp_shell = 0
        if( cline%defined('lp') ) lp_shell = calc_fourier_index(params%lp, params%box, params%smpd)
        hp_shell = calc_fourier_index(params%hp, params%box, params%smpd)
        call series%set_band(hp_shell, lp_shell)
        call series%fit
        call series%get_band(klo, khi)
        write(logfhandle,'(A,I0,A,I0)') '>>> DOSE SERIES: FITTING BAND SHELLS ', klo, ' TO ', khi
        call series%paired_cv(em, om, shell)
        deallocate(em, om)
        ! pass two: the halves as given, accumulated under the model, the control and the flat mean,
        ! beside window 1 alone, so the gain is on one line
        call gather_series(win, params, ncoeff, .false., em, om)
        allocate(ae(ncoeff), ao(ncoeff), fsc(nsh), fsc_msk(nsh), fsc_tab(nsh,4))
        outbody = strip_extension(params%outvol)
        res     = vol%get_res()
        call series%accumulate(em, om, shell, DOSE_WEIGHTS_MODEL, ae, ao)
        call finish_pair(ae, ao, params, fsc, fsc_msk, outbody, .true.)
        call get_resolution(fsc,     res, fsc05(1),  fsc0143(1))
        call get_resolution(fsc_msk, res, mfsc05(1), mfsc0143(1))
        fsc_tab(:,1) = fsc
        fsc_tab(:,2) = fsc_msk
        if( params%zerodose .eq. 'yes' )then
            allocate(s0(ncoeff))
            call series%zerodose(0.5*(ae + ao), shell, s0)
            call write_coefficients(s0, params, outbody//'_zerodose.mrc')
            deallocate(s0)
        endif
        call series%accumulate(em, om, shell, DOSE_WEIGHTS_CONTROL, ae, ao)
        call finish_pair(ae, ao, params, fsc, fsc_msk, outbody//'_control', .true.)
        call get_resolution(fsc,     res, fsc05(2),  fsc0143(2))
        call get_resolution(fsc_msk, res, mfsc05(2), mfsc0143(2))
        fsc_tab(:,3) = fsc
        fsc_tab(:,4) = fsc_msk
        call series%accumulate(em, om, shell, DOSE_WEIGHTS_FLAT, ae, ao)
        call finish_pair(ae, ao, params, fsc, fsc_msk, outbody, .false.)
        call get_resolution(fsc,     res, fsc05(3),  fsc0143(3))
        call get_resolution(fsc_msk, res, mfsc05(3), mfsc0143(3))
        call finish_pair(em(:,1), om(:,1), params, fsc, fsc_msk, outbody, .false.)
        call get_resolution(fsc,     res, fsc05(4),  fsc0143(4))
        call get_resolution(fsc_msk, res, mfsc05(4), mfsc0143(4))
        deallocate(em, om, ae, ao)
        write(logfhandle,'(A,4(A,2F8.2))') '>>> Dose series: FSC=0.5 and 0.143 in A', &
            &'  model', fsc05(1), fsc0143(1), '  control', fsc05(2), fsc0143(2), &
            &'  flat', fsc05(3), fsc0143(3), '  window 1', fsc05(4), fsc0143(4)
        write(logfhandle,'(A,4(A,2F8.2))') '>>> Dose series: masked FSC=0.5 and 0.143 in A', &
            &'  model', mfsc05(1), mfsc0143(1), '  control', mfsc05(2), mfsc0143(2), &
            &'  flat', mfsc05(3), mfsc0143(3), '  window 1', mfsc05(4), mfsc0143(4)
        call write_fsc_table(outbody//'_fsc.txt', res, fsc_tab)
        call series%write_model_table(outbody//'_model.txt')
        call series%write_attenuation_table(outbody//'_attenuation.txt')
        call series%write_weights_table(outbody//'_weights.txt')
        call series%write_cv_table(outbody//'_cv.txt')
        call series%write_decay_table(outbody//'_decay.txt')
        call series%write_step_exposure(outbody//'_decay_steps.txt', PLOT_RES)
        call series%write_resolution_weights(outbody//'_weights_res.txt', PLOT_RES)
        ! the model as reextract reads it; the manifest's windows share one physical dose per fraction
        call series%export_curve(curve, win(:)%ffrom, win(:)%fto, &
            &(win(1)%dose_hi - win(1)%dose_lo) / real(win(1)%fto - win(1)%ffrom + 1, dp))
        call curve%write(outbody//'_curve.txt')
        call plot_decay(outbody)
        call series%log_readout
        ! destruct
        call curve%kill
        call series%kill
        call vol%kill
        call mskimg%kill
        call simple_end('**** SIMPLE_DOSE_FILTER3D NORMAL STOP ****')

        contains

            !> the output body: outvol up to the extension of its file name; a dot in a directory is kept
            function strip_extension( fname ) result( body )
                class(string), intent(in) :: fname
                type(string) :: body
                character(len=:), allocatable :: txt
                integer :: pos
                txt = fname%to_char()
                pos = index(txt, '.', back=.true.)
                if( pos > 1 .and. pos > index(txt, '/', back=.true.) + 1 )then
                    body = txt(1:pos-1)
                else
                    body = txt
                endif
            end function strip_extension

            !> shell index of every coefficient of the half-volume lattice, in a fixed traversal order
            subroutine lattice_shells( vol, nsh, shell, ncoeff )
                class(image),         intent(in)    :: vol
                integer,              intent(in)    :: nsh
                integer, allocatable, intent(inout) :: shell(:)
                integer,              intent(out)   :: ncoeff
                integer :: lims(3,2), h, k, l, s, i
                lims   = vol%loop_lims(2)
                ncoeff = 0
                do l = lims(3,1), lims(3,2)
                    do k = lims(2,1), lims(2,2)
                        do h = lims(1,1), lims(1,2)
                            ncoeff = ncoeff + 1
                        enddo
                    enddo
                enddo
                if( allocated(shell) ) deallocate(shell)
                allocate(shell(ncoeff))
                i = 0
                do l = lims(3,1), lims(3,2)
                    do k = lims(2,1), lims(2,2)
                        do h = lims(1,1), lims(1,2)
                            i = i + 1
                            ! DC is 0 and the cube's corners past Nyquist nsh+1: both stay out of the statistics,
                            ! and the accumulation clamps them to the first and last shell's weights
                            s = nint(sqrt(real(h*h + k*k + l*l)))
                            shell(i) = min(s, nsh + 1)
                        enddo
                    enddo
                enddo
            end subroutine lattice_shells

            !> read every window's halves and gather their coefficients; l_mask selects the statistics copy
            subroutine gather_series( win, params, ncoeff, l_mask, em, om )
                type(dose_window),        intent(in)    :: win(:)
                class(parameters),        intent(in)    :: params
                integer,                  intent(in)    :: ncoeff
                logical,                  intent(in)    :: l_mask
                complex(sp), allocatable, intent(inout) :: em(:,:), om(:,:)
                type(image) :: vol
                integer     :: w, nwin
                nwin = size(win)
                if( allocated(em) ) deallocate(em)
                if( allocated(om) ) deallocate(om)
                allocate(em(ncoeff,nwin), om(ncoeff,nwin), source=cmplx(0.,0.))
                call vol%new([params%box,params%box,params%box], params%smpd)
                do w = 1, nwin
                    call prepare_half(vol, win(w)%even, l_mask, params)
                    call gather_half(vol, em(:,w))
                    call prepare_half(vol, win(w)%odd, l_mask, params)
                    call gather_half(vol, om(:,w))
                enddo
                call vol%kill
            end subroutine gather_series

            !> read a half and transform it: soft-masked about its solvent level for the statistics,
            !! mean removed and otherwise as given for the accumulation
            subroutine prepare_half( vol, fname, l_mask, params )
                class(image),      intent(inout) :: vol
                class(string),     intent(in)    :: fname
                logical,           intent(in)    :: l_mask
                class(parameters), intent(in)    :: params
                ! Headers were checked against the native grid before any image buffer was allocated.
                ! the image is reused across reads: a Fourier flag left by the last fft would skip the next one
                call vol%zero_and_unflag_ft
                call vol%read(fname)
                if( l_mask )then
                    call vol%mask3D_soft(params%msk)
                else
                    call vol%subtr(vol%mean())
                endif
                call vol%fft()
            end subroutine prepare_half

            !> pull the half-volume lattice into a flat array in the order lattice_shells fixed
            subroutine gather_half( vol, dest )
                class(image), intent(in)  :: vol
                complex(sp),  intent(out) :: dest(:)
                integer :: lims(3,2), phys(3), h, k, l, i
                lims = vol%loop_lims(2)
                i    = 0
                do l = lims(3,1), lims(3,2)
                    do k = lims(2,1), lims(2,2)
                        do h = lims(1,1), lims(1,2)
                            i    = i + 1
                            phys = vol%comp_addr_phys(h,k,l)
                            dest(i) = vol%get_cmat_at(phys(1),phys(2),phys(3))
                        enddo
                    enddo
                enddo
            end subroutine gather_half

            !> write a flat coefficient array back onto the lattice; the volume is left in Fourier space
            subroutine scatter_half( vol, src )
                class(image), intent(inout) :: vol
                complex(sp),  intent(in)    :: src(:)
                integer :: lims(3,2), phys(3), h, k, l, i
                call vol%zero_and_flag_ft
                lims = vol%loop_lims(2)
                i    = 0
                do l = lims(3,1), lims(3,2)
                    do k = lims(2,1), lims(2,2)
                        do h = lims(1,1), lims(1,2)
                            i    = i + 1
                            phys = vol%comp_addr_phys(h,k,l)
                            call vol%set_cmat_at(phys(1),phys(2),phys(3), src(i))
                        enddo
                    enddo
                enddo
            end subroutine scatter_half

            !> scatter an accumulated pair and take its FSC in Fourier space; optionally write the halves,
            !! their FSC and their mean under body
            !> the FSC of a half pair, unmasked and under the soft mask of mskdiam; with l_write the halves,
            !! their mean and both FSCs are written
            subroutine finish_pair( ae, ao, params, fsc, fsc_msk, body, l_write )
                complex(sp),       intent(in)    :: ae(:), ao(:)
                class(parameters), intent(in)    :: params
                real,              intent(inout) :: fsc(:), fsc_msk(:)
                class(string),     intent(in)    :: body
                logical,           intent(in)    :: l_write
                type(image) :: he, ho, me, mo
                call he%new([params%box,params%box,params%box], params%smpd)
                call ho%new([params%box,params%box,params%box], params%smpd)
                call scatter_half(he, ae)
                call scatter_half(ho, ao)
                call he%fsc(ho, fsc)
                call he%ifft()
                call ho%ifft()
                call me%copy(he)
                call mo%copy(ho)
                call me%mask3D_soft(params%msk, backgr=0.)
                call mo%mask3D_soft(params%msk, backgr=0.)
                call me%fft()
                call mo%fft()
                call me%fsc(mo, fsc_msk)
                if( l_write )then
                    call he%write(body//'_even.mrc')
                    call ho%write(body//'_odd.mrc')
                    call arr2file(fsc, body//'_fsc.bin')
                    call arr2file(fsc_msk, body//'_fsc_masked.bin')
                    call he%add(ho)
                    call he%mul(0.5)
                    call he%write(body//'.mrc')
                endif
                call he%kill
                call ho%kill
                call me%kill
                call mo%kill
            end subroutine finish_pair

            !> per shell: resolution, then the model's and the control's FSC, each unmasked and masked
            subroutine write_fsc_table( fname, res, tab )
                class(string), intent(in) :: fname
                real,          intent(in) :: res(:), tab(:,:)
                character, parameter :: TAB_CHAR = achar(9)
                integer :: funit, s, io_stat
                call fopen(funit, file=fname, status='replace', action='write', form='formatted')
                write(funit,'(A)',iostat=io_stat) '# FSC of the even and odd accumulations; masked under the soft mask of mskdiam'
                write(funit,'(A)',iostat=io_stat) 'shell'//TAB_CHAR//'res_A'//TAB_CHAR//'model'//TAB_CHAR//'model_masked'// &
                    &TAB_CHAR//'control'//TAB_CHAR//'control_masked'
                do s = 1, size(res)
                    write(funit,'(I0,5(A,F10.5))',iostat=io_stat) s, TAB_CHAR, res(s), TAB_CHAR, tab(s,1), TAB_CHAR, tab(s,2), &
                        &TAB_CHAR, tab(s,3), TAB_CHAR, tab(s,4)
                enddo
                call fileiochk('Writing the FSC table '//fname%to_char(), io_stat)
                call fclose(funit)
            end subroutine write_fsc_table

            !> scatter a flat coefficient array and write it as a real-space map
            subroutine write_coefficients( coeffs, params, fname )
                complex(sp),       intent(in) :: coeffs(:)
                class(parameters), intent(in) :: params
                class(string),     intent(in) :: fname
                type(image) :: map
                call map%new([params%box,params%box,params%box], params%smpd)
                call scatter_half(map, coeffs)
                call map%ifft()
                call map%write(fname)
                call map%kill
            end subroutine write_coefficients

            !> The decay figure (fading per window, critical exposure per shell and between neighbouring windows,
            !! all against Grant and Grigorieff) and the weights at PLOT_RES, rendered by gnuplot from the tables
            !! beside them; a failure warns and keeps the tables.
            subroutine plot_decay( body )
                class(string), intent(in) :: body
                character(len=:), allocatable :: b, gpfile, logfile, line, msg
                integer :: unit, ios, i, j, s, w, ib, nsel, sel(5), cmdstat, exitstat
                logical :: l_decay, l_weights
                b       = body%to_char()
                gpfile  = b//'_plots.gp'
                logfile = b//'_plots_gnuplot.log'
                ! up to five windows after the first, evenly spaced, the last among them
                nsel = min(5, nwin-1)
                do j = 1, nsel
                    sel(j) = 2 + nint(real((j-1)*(nwin-2),dp) / real(max(1,nsel-1),dp))
                enddo
                open(newunit=unit, file=gpfile, status='replace', action='write', iostat=ios)
                if( ios /= 0 )then
                    THROW_WARN('Cannot write the plot commands; the tables are kept')
                    return
                endif
                call gp(unit, 'set terminal pngcairo size 2600,800 enhanced font "Sans,11"')
                call gp(unit, 'set datafile missing "NaN"')
                call gp(unit, 'set output '//gq(b//'_decay.png'))
                call gp(unit, 'set multiplot layout 1,3 title "Dose series: measured fading against Grant and Grigorieff"')
                call gp(unit, 'set title "Attenuation relative to window 1: points measured, solid model, dashed Grant and Grigorieff"')
                call gp(unit, 'set xlabel "Spatial frequency (Angstrom^{-1})"; set ylabel "a_w(s)"; set yrange [-0.1:1.2]')
                call gp(unit, 'set key top right font ",9"')
                line = 'plot '
                do j = 1, nsel
                    w  = sel(j)
                    ib = 7 + 4*(w-1)
                    if( j > 1 ) line = line//', '
                    line = line//gq(b//'_decay.txt')//' using 2:'//int2str(ib+1)//':'//int2str(ib+2)// &
                        &' with yerrorbars pt 7 ps 0.6 lc '//int2str(j)//' title '// &
                        &gq('fractions '//int2str(win(w)%ffrom)//'-'//int2str(win(w)%fto))// &
                        &', '''' using 2:'//int2str(ib+3)//' with lines lw 2 lc '//int2str(j)//' notitle'// &
                        &', '''' using 2:'//int2str(ib+4)//' with lines dt 2 lw 1.5 lc '//int2str(j)//' notitle'
                enddo
                call gp(unit, line)
                call gp(unit, 'set title "Critical exposure per shell"; set logscale y; set yrange [*:*]')
                call gp(unit, 'set ylabel "N_e (e/Angstrom^2, scaled dose)"; set key top right font ",10"')
                call gp(unit, 'plot '//gq(b//'_decay.txt')//' using 2:3:4 with yerrorbars pt 7 ps 0.6 lc 1 '// &
                    &'title "measured per shell, no model", '''' using 2:7 with lines lw 2 lc 2 '// &
                    &'title "model per window, same estimator", '''' using 2:5 with lines dt 3 lw 2 lc 2 '// &
                    &'title "model lines, 1/(2 gamma + beta k^2/2)", '''' using 2:6 with lines dt 2 lw 2 lc 3 '// &
                    &'title "Grant and Grigorieff"')
                ! flat lines here would be fading exponential in dose, as the lines and Grant and Grigorieff assume
                call gp(unit, 'set title "Critical exposure between neighbouring windows: points measured, '// &
                    &'dashed Grant and Grigorieff"')
                call gp(unit, 'set xlabel "Mean dose of the two windows (e/Angstrom^2, scaled)"; set key top right font ",10"')
                line = 'plot '
                j = 0
                do i = 1, size(PLOT_RES)
                    s = nint(1._dp / (PLOT_RES(i) * kfreq(1)))
                    if( s < klo .or. s > khi ) cycle    ! outside the fitting band
                    j = j + 1
                    if( j > 1 ) line = line//', '
                    line = line//gq(b//'_decay_steps.txt')//' using 2:'//int2str(3*i)//':'//int2str(3*i+1)// &
                        &' with yerrorlines pt 7 lw 2 lc '//int2str(j)//' title '//gq(int2str(nint(PLOT_RES(i)))//' A')// &
                        &', '''' using 2:'//int2str(3*i+2)//' with lines dt 2 lw 1.5 lc '//int2str(j)//' notitle'
                enddo
                if( j > 0 ) call gp(unit, line)
                call gp(unit, 'unset multiplot; unset logscale y')
                call gp(unit, 'set terminal pngcairo size 1100,800 enhanced font "Sans,11"')
                call gp(unit, 'set output '//gq(b//'_weights_res.png'))
                call gp(unit, 'set title "Accumulation weights per window: solid model, dashed Grant and Grigorieff"')
                call gp(unit, 'set xlabel "Window dose centre (e/Angstrom^2, scaled)"; set ylabel "weight"; set yrange [0:*]')
                call gp(unit, 'set key top right font ",10"')
                line = 'plot '
                j = 0
                do i = 1, size(PLOT_RES)
                    if( PLOT_RES(i) < 2._dp * real(params%smpd,dp) ) cycle    ! past Nyquist
                    j = j + 1
                    if( j > 1 ) line = line//', '
                    line = line//gq(b//'_weights_res.txt')//' using 2:'//int2str(1+2*i)// &
                        &' with linespoints lw 2 pt 7 lc '//int2str(j)//' title '//gq(int2str(nint(PLOT_RES(i)))//' A')// &
                        &', '''' using 2:'//int2str(2+2*i)//' with lines dt 2 lw 1.5 lc '//int2str(j)//' notitle'
                enddo
                line = line//', 1./'//int2str(nwin)//' with lines dt 3 lc "black" title "equal weights"'
                call gp(unit, line)
                call gp(unit, 'unset output')
                close(unit)
                if( file_exists(b//'_decay.png')       ) call del_file(b//'_decay.png')
                if( file_exists(b//'_weights_res.png') ) call del_file(b//'_weights_res.png')
                exitstat = 1; cmdstat = 0
                call execute_command_line('gnuplot '//posix_quote(gpfile)//' > '//posix_quote(logfile)//' 2>&1', &
                    &wait=.true., exitstat=exitstat, cmdstat=cmdstat)
                l_decay   = file_exists(b//'_decay.png')
                l_weights = file_exists(b//'_weights_res.png')
                if( cmdstat /= 0 .or. exitstat /= 0 .or. .not. l_decay .or. .not. l_weights )then
                    msg = 'Dose series plots not rendered; the tables are kept, see '//logfile
                    THROW_WARN(msg)
                else
                    write(logfhandle,'(A)') '>>> DOSE SERIES: WROTE '//b//'_decay.png AND '//b//'_weights_res.png'
                endif
            end subroutine plot_decay

            !> one gnuplot command line
            subroutine gp( unit, command )
                integer,          intent(in) :: unit
                character(len=*), intent(in) :: command
                integer :: io
                write(unit,'(A)',iostat=io) command
            end subroutine gp

            !> text as a gnuplot single-quoted string
            function gq( text ) result( quoted )
                character(len=*), intent(in)  :: text
                character(len=:), allocatable :: quoted
                quoted = trim(fortran_quote(text))
            end function gq

    end subroutine exec_dose_filter3D

    !> calculates Fourier shell correlation from Even/Odd Volume pairs
    subroutine exec_fsc( self, cline )
        class(commander_fsc), intent(inout) :: self
        class(cmdline),       intent(inout) :: cline
        real, allocatable :: res(:), fsc(:), fsc_t(:), fsc_n(:)
        type(parameters)  :: params
        type(image)       :: even, odd, merged
        type(image_msk)   :: mskvol
        type(string)      :: fsc_templ
        integer           :: j, k_hp, k_lp, nyq
        real              :: res_fsc05, res_fsc0143
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        ! read even/odd pair
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        res = even%get_res()
        nyq = even%get_filtsz()
        if( params%automsk .ne. 'no' )then
            call merged%copy(even)
            call merged%add(odd)
            call merged%mul(0.5)
            call mskvol%automask3D(params, merged, params%automsk.eq.'tight')
            call phase_rand_fsc(even, odd, mskvol, 1, nyq, fsc, fsc_t, fsc_n)
            call merged%kill
        else
            ! spherical masking
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
            call even%fft()
            call odd%fft()
            allocate(fsc(nyq),source=0.)
            call even%fsc(odd, fsc)
        endif
        do j=1,nyq
            write(logfhandle,'(A,1X,F6.2,1X,A,1X,F7.3)') '>>> RESOLUTION:', res(j), '>>> CORRELATION:', fsc(j)
        end do
        call get_resolution(fsc, res, res_fsc05, res_fsc0143)
        k_hp = calc_fourier_index(params%hp, params%box, params%smpd)
        k_lp = calc_fourier_index(params%lp, params%box, params%smpd)
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.500 DETERMINED TO:', res_fsc05
        write(logfhandle,'(A,1X,F6.2)') '>>> RESOLUTION AT FSC=0.143 DETERMINED TO:', res_fsc0143
        call even%kill
        call odd%kill
        fsc_templ = refine3D_fsc_fbody(1)
        call arr2file(fsc, refine3D_fsc_fname(1))
        if( params%automsk .ne. 'no' )then
            call plot_phrand_fsc(size(fsc), fsc, fsc_t, fsc_n, res, params%smpd, fsc_templ%to_char())
        else
            call plot_fsc(size(fsc), fsc, res, params%smpd, fsc_templ%to_char())
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_FSC NORMAL STOP ****')
    end subroutine exec_fsc

    !> calculates a CryoSPARC-like conical FSC area ratio from Even/Odd Volume pairs
    subroutine exec_fsc_area_score( self, cline )
        class(commander_fsc_area_score), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)             :: params
        type(image)                  :: even, odd, merged
        type(image_msk)              :: mskvol
        type(fsc_area_score_result)  :: result
        type(string)                 :: fbody
        if( .not. cline%defined('mkdir')      ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('athres')     ) call cline%set('athres', 20.)
        if( .not. cline%defined('nspace')     ) call cline%set('nspace', 256)
        if( .not. cline%defined('lplim_crit') ) call cline%set('lplim_crit', 0.143)
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        if( params%automsk .ne. 'no' )then
            call merged%copy(even)
            call merged%add(odd)
            call merged%mul(0.5)
            call mskvol%automask3D(params, merged, params%automsk.eq.'tight')
            call mskvol%apply_3Dmask(even)
            call mskvol%apply_3Dmask(odd)
        else
            call even%mask3D_soft(params%msk)
            call odd%mask3D_soft(params%msk)
        endif
        call even%fft()
        call odd%fft()
        call result%new(even, params%nspace, params%athres, params%lplim_crit, 1)
        call result%calc_fsc_area_score(even, odd)
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = 'fsc_area_score'
        endif
        call write_fsc_area_score_outputs(result, fbody%to_char())
        call plot_fsc_area_score(result, fbody%to_char())
        call even%kill
        call odd%kill
        call mskvol%kill_bimg
        call result%kill
        call merged%kill
        call simple_end('**** SIMPLE_FSC_AREA_SCORE NORMAL STOP ****')
    end subroutine exec_fsc_area_score

    subroutine exec_uniform_filter3D(self, cline)
        use simple_opt_filter, only: estimate_lplim
        class(commander_uniform_filter3D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(parameters)  :: params
        type(image)       :: merged, even, odd, odd_filt
        type(image_msk)   :: mskvol
        real              :: lpopt
        type(string)      :: file_tag
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call odd_filt%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd %read(params%vols(1))
        call even%read(params%vols(2))
        file_tag = 'uniform_3D_filter'
        if( params%automsk .ne. 'no' )then
            call merged%copy(even)
            call merged%add(odd)
            call merged%mul(0.5)
            call mskvol%automask3D(params, merged, l_tight=params%automsk.eq.'tight')
            call mskvol%one_at_edge ! to expand before masking of reference internally
            call merged%kill
        else
            call mskvol%disc([params%box,params%box,params%box], params%smpd,&
                    &real(min(params%box/2, int(params%msk + COSMSKHALFWIDTH))))
        endif
        ! soft masking needed for FT
        call even%mask3D_soft(params%msk)
        call  odd%mask3D_soft(params%msk)
        call estimate_lplim(odd, even, mskvol, params%lpstart, params%lpstop, lpopt, odd_filt)
        print *, 'found optimal low-pass limit: ', lpopt
        call odd_filt%write(string('odd_filt.mrc'))
        ! destruct
        call odd%kill
        call odd_filt%kill
        call even%kill
        call mskvol%kill_bimg
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER3D NORMAL STOP ****')
    end subroutine exec_uniform_filter3D

    subroutine exec_nu_filt3D(self, cline)
        use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter, &
            &print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, write_nu_local_resolution_map, NU_DEV_OUTPUT,&
            &nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, write_nu_evidence_map, print_nu_envmask_stats,&
            &NU_ENVMASK_MINVOL_FRAC, NU_ENVMASK_GROW_A, NU_ENVMASK_EDGE_A
        class(commander_nu_filt3D), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)     :: params
        type(image)          :: even, odd, even_nu, odd_nu
        type(image_msk)      :: nu_envelope
        type(nu_envmask_params) :: envp
        type(nu_envmask_stats)  :: envstats
        type(string)         :: even_out, odd_out, avg_out, locres_out, evid_out, envmsk_out
        logical, allocatable :: l_env(:,:,:)
        integer              :: n_ccs, n_ccs_kept, grow_px, edge_px
        logical              :: l_envmask
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call odd%new([params%box,params%box,params%box], params%smpd)
        call even%new([params%box,params%box,params%box], params%smpd)
        call odd%read(params%vols(1))
        call even%read(params%vols(2))
        call setup_nu_dmats(even, odd, params%mskdiam, [real ::])
        call optimize_nu_cutoff_finds()
        odd_out  = add2fbody(params%vols(1), params%ext, NUFILT_SUFFIX)
        even_out = add2fbody(params%vols(2), params%ext, NUFILT_SUFFIX)
        if( params%outvol .ne. '' )then
            avg_out = params%outvol
        else
            avg_out = 'vol_nufilt'//params%ext%to_char()
        endif
        locres_out = add2fbody(avg_out, params%ext, NULOCRES_SUFFIX)
        ! The envelope reads the raw evidence (nu_ev_base, nu_ev_best), which
        ! lives until cleanup_nu_filter.
        l_envmask = params%nu_envmsk .eq. 'yes'
        if( l_envmask )then
            ! The evidence margin is always written: its histogram is what tells
            ! you whether the solvent null separates from the density mode at all,
            ! which is the precondition for the mask being meaningful.
            envp%nsigma      = params%nu_msk_sig
            envp%lp_smooth   = params%amsklp
            ! beta, dens_weight and l_relative stay at the NU_ENVMASK_* constants of
            ! the type definition, as in refinement.
            evid_out = add2fbody(avg_out, params%ext, NUEVIDENCE_SUFFIX)
            call write_nu_evidence_map(evid_out, envp%lp_smooth, envp%l_relative)
            call nu_evidence_envelope(envp, l_env, envstats)
            call print_nu_envmask_stats(envstats)
            grow_px = max(1, nint(NU_ENVMASK_GROW_A / even%get_smpd()))
            edge_px = max(1, nint(NU_ENVMASK_EDGE_A / even%get_smpd()))
            call nu_envelope%envmask3D_from_lmask(l_env, even%get_smpd(), grow_px, edge_px, &
                &NU_ENVMASK_MINVOL_FRAC, .true., n_ccs, n_ccs_kept)
            ! an empty evidence field returns without constructing the mask
            ! image; writing it would abort on invalid MRC dimensions
            if( n_ccs_kept < 1 )then
                THROW_WARN('NU evidence envelope is empty at these settings; no envelope mask written. Feed the _unfil half pair when the inputs come from an ml_reg run (regularized maps flatten the evidence margin), or loosen nu_msk_sig / amsklp')
            else
                envmsk_out = add2fbody(avg_out, params%ext, NUENVMSK_SUFFIX)
                call nu_envelope%write(envmsk_out, del_if_exists=.true.)
                call wait_for_closure(envmsk_out)
            endif
        endif
        call nu_filter_vols(even_nu, odd_nu, even, odd)
        call print_nu_filtmap_lowpass_stats()
        if( NU_DEV_OUTPUT ) call analyze_filtmap_neighbor_continuity()
        call odd_nu%write(odd_out, del_if_exists=.true.)
        call even_nu%write(even_out, del_if_exists=.true.)
        call even_nu%add(odd_nu)
        call even_nu%mul(0.5)
        call even_nu%write(avg_out, del_if_exists=.true.)
        call write_nu_local_resolution_map(locres_out)
        call wait_for_closure(avg_out)
        call wait_for_closure(locres_out)
        call cleanup_nu_filter()
        call odd_nu%kill
        call even_nu%kill
        call odd%kill
        call even%kill
        call nu_envelope%kill_bimg
        call locres_out%kill
        call evid_out%kill
        call envmsk_out%kill
        if( allocated(l_env)  ) deallocate(l_env)
        call simple_end('**** SIMPLE_nu_filt3D NORMAL STOP ****')
    end subroutine exec_nu_filt3D

    subroutine exec_uniform_filter2D( self, cline )
        use simple_opt_filter, only: estimate_lplims2D
        class(commander_uniform_filter2D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(image),      allocatable :: odd(:), even(:), odd_filt(:)
        real,             allocatable :: lpsopt(:)
        type(parameters) :: params
        integer          :: iptcl
        ! init
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call params%new(cline)
        call find_ldim_nptcls(params%stk, params%ldim, params%nptcls)
        params%ldim(3) = 1 ! because we operate on stacks
        ! allocate
        allocate(odd(params%nptcls), odd_filt(params%nptcls), even(params%nptcls), lpsopt(params%nptcls))
        ! construct & read
        do iptcl = 1, params%nptcls
            call odd( iptcl)%new( params%ldim, params%smpd, .false.)
            call odd_filt( iptcl)%new( params%ldim, params%smpd, .false.)
            call odd( iptcl)%read(params%stk,  iptcl)
            call even(iptcl)%new( params%ldim, params%smpd, .false.)
            call even(iptcl)%read(params%stk2, iptcl)
        enddo
        ! filter
        call estimate_lplims2D( odd, even, params%msk, params%lpstart, params%lpstop, lpsopt, odd_filt )
        ! write output and destruct
        do iptcl = 1, params%nptcls
            print *, 'found optimal low-pass limit: ', iptcl, lpsopt(iptcl)
            call odd_filt(iptcl)%write(string('odd_filt.mrc'), iptcl)
        end do
        ! end gracefully
        call simple_end('**** SIMPLE_UNIFORM_FILTER2D NORMAL STOP ****')
    end subroutine exec_uniform_filter2D

    subroutine exec_estimate_lpstages( self, cline )
        class(commander_estimate_lpstages), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(builder)    :: build
        type(parameters) :: params
        real, parameter  :: LPSTART_LB=10., LPSTART_DEFAULT=20.
        real, parameter  :: LPSTOP_BOUNDS(2)  = [4.5,6.0]
        real, parameter  :: LPSTART_BOUNDS(2) = [10.,20.] 
        type(string)     :: frcs_fname
        real,              allocatable :: frcs_avg(:)
        integer,           allocatable :: states(:)
        type(lp_crop_inf), allocatable :: lpinfo(:)
        integer :: filtsz
        real    :: lpfinal
        call build%init_params_and_build_general_tbox(cline, params)
        call build%spproj%get_frcs(frcs_fname, 'frc2D', fail=.true.)
        call build%clsfrcs%read(frcs_fname)
        states  = nint(build%spproj%os_cls2D%get_all('state'))
        filtsz = build%clsfrcs%get_filtsz()
        allocate(frcs_avg(filtsz), source=0.)
        call build%clsfrcs%avg_frc_getter(frcs_avg, states)
        allocate(lpinfo(params%nstages))
        lpfinal = max(LPSTOP_BOUNDS(1),calc_lplim_final_stage(3))
        lpfinal = min(LPSTOP_BOUNDS(2),lpfinal)
        call lpstages(params%box, params%nstages, frcs_avg, params%smpd,&
        &LPSTART_BOUNDS(1), LPSTART_BOUNDS(2), lpfinal, lpinfo, l_cavgs=.false.)
        call simple_end('**** SIMPLE_ESTIMATE_LPSTAGES NORMAL STOP ****')

        contains

            function calc_lplim_final_stage( nbest ) result( lplim )
                integer, intent(in)  :: nbest
                real,    allocatable :: res(:), tmp_rarr(:)
                integer, allocatable :: tmp_iarr(:)
                real :: lplim
                tmp_rarr  = build%spproj%os_cls2D%get_all('res')
                tmp_iarr  = nint(build%spproj%os_cls2D%get_all('state'))
                res       = pack(tmp_rarr, mask=(tmp_iarr>0))
                call hpsort(res)
                lplim = median_nocopy(res(:nbest))
                deallocate(tmp_rarr, tmp_iarr, res)
            end function calc_lplim_final_stage

    end subroutine exec_estimate_lpstages

end module simple_commanders_resolest
