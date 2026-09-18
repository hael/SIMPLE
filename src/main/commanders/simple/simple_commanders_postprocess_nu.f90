!@descr: NU-evidence nonuniform postprocessing (isolated from the standard postprocess path)
!
! postprocess_nu is the commander for the NU-evidence local sharpening
! experiment (nu_evidence_local_sharpening.md): model-free LocScale-style
! local amplitude restoration in which both the confidence field and the
! target spectrum derive from the frozen cross-half NU evidence state (the
! compact state the NU competition also uses for its envelope). It is deliberately
! isolated from the standard postprocess commander (global B-factor + FSC
! filter), which remains untouched: a single isotropic B-factor does not
! serve most specimens, and this path is the recorded alternative.
!
! Operates on the project like postprocess (2026-09-18): the state's volume
! from the out segment, its UNREGULARIZED even/odd pair (_even_unfil/_odd_unfil;
! evidence authority is the base pair, regularized maps flatten the evidence
! margin) as the evidence input, and its regularized pair (_even/_odd) as the
! auxiliary member of the refinement's filter competition, whose products are
! written first. Every product carries the _pproc_nu suffix (the sharpened
! map: _pproc_nu; the competition's references and local-resolution map:
! _pproc_nu_filt, _pproc_nu_locres) and is a display/interpretation map --
! never an input to FSC correction or resolution claims.
module simple_commanders_postprocess_nu
use simple_commanders_api
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, &
    &get_nu_filter_bank_finest_lp, build_nu_evidence_state, cleanup_nu_filter, nu_evidence_state, &
    &assert_nu_evidence_replay_ready, print_nu_evidence_summary, nu_evidence_sharpen_vol, NU_EVIDENCE_SOURCE_BASE, &
    &nu_filter_vols, print_nu_filtmap_lowpass_stats, write_nu_local_resolution_map, &
    &get_nu_filtmap_finest_selected_lp, NU_ALIGN_LP_MIN_SIGNAL_PCT
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_postprocess_nu
  contains
    procedure :: execute      => exec_postprocess_nu
end type commander_postprocess_nu

contains

    subroutine exec_postprocess_nu( self, cline )
        class(commander_postprocess_nu), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)        :: params
        type(sp_project)        :: spproj
        type(nu_evidence_state) :: evstate
        type(image)             :: even, odd, vol_sharp, vol_even_nu, vol_odd_nu
        type(image), allocatable :: aux_even(:), aux_odd(:)
        type(string)            :: vol_out, fname, fname_vol, fname_even_unfil, fname_odd_unfil, fname_even, fname_odd
        real, allocatable       :: corrs(:), res(:)
        real                    :: fsc05, fsc0143, handoff_lp, raw_lp, smpd
        integer                 :: n_signal, state, box, ldim(3), nptcls
        logical                 :: l_aux
        ! operates on the project like postprocess: the state's volume from
        ! the out segment, its unregularized pair (_even_unfil/_odd_unfil)
        ! as the evidence input and its regularized pair (_even/_odd) as the
        ! auxiliary member of the filter competition when present
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
        call cline%set('oritype', 'out')
        call cline%delete('ptcl_src')
        call params%new(cline)
        call spproj%read_segment(params%oritype, params%projfile)
        state = 1
        if( cline%defined('state') ) state = params%state
        call spproj%get_vol('vol', state, fname_vol, smpd, box)
        call spproj%kill
        if( .not. file_exists(fname_vol) ) THROW_HARD('volume: '//fname_vol%to_char()//' does not exist')
        fname_even_unfil = add2fbody(fname_vol, params%ext, '_even_unfil')
        fname_odd_unfil  = add2fbody(fname_vol, params%ext, '_odd_unfil')
        fname_even       = add2fbody(fname_vol, params%ext, '_even')
        fname_odd        = add2fbody(fname_vol, params%ext, '_odd')
        if( .not. file_exists(fname_even_unfil) .or. .not. file_exists(fname_odd_unfil) ) &
            &THROW_HARD('the unregularized half pair (_even_unfil/_odd_unfil) of '//fname_vol%to_char()//' is missing; postprocess_nu')
        call find_ldim_nptcls(fname_even_unfil, ldim, nptcls)
        if( ldim(1) /= box ) THROW_HARD('the unregularized pair is at a different box than the project volume; postprocess_nu')
        params%box  = box
        params%smpd = smpd
        write(logfhandle,'(A,I0,A,I0,A,F7.3)') '>>> POSTPROCESS_NU: STATE ', state, ', box ', box, ', smpd ', smpd
        write(logfhandle,'(A)') '>>> POSTPROCESS_NU: UNREGULARIZED PAIR '//fname_odd_unfil%to_char()//' '//fname_even_unfil%to_char()
        call odd%new([box,box,box], smpd)
        call even%new([box,box,box], smpd)
        call odd%read(fname_odd_unfil)
        call even%read(fname_even_unfil)
        l_aux = file_exists(fname_even) .and. file_exists(fname_odd)
        if( l_aux )then
            call find_ldim_nptcls(fname_even, ldim, nptcls)
            l_aux = ldim(1) == box
        endif
        if( l_aux )then
            write(logfhandle,'(A)') '>>> POSTPROCESS_NU: REGULARIZED PAIR '//fname_odd%to_char()//' '//fname_even%to_char()
            allocate(aux_odd(1), aux_even(1))
            call aux_odd(1)%new([box,box,box], smpd)
            call aux_even(1)%new([box,box,box], smpd)
            call aux_odd(1)%read(fname_odd)
            call aux_even(1)%read(fname_even)
        else
            write(logfhandle,'(A)') '>>> POSTPROCESS_NU: no regularized pair beside the project volume; filter competition skipped'
        endif
        ! the same bank as the refinement: the static ladder capped at the
        ! pair's FSC=0.143 / NU_BANK_FSC_HEADROOM; the evidence compaction
        ! takes the base pair's candidates only (its contract)
        allocate(corrs(fdim(box)-1), source=0.)
        ! image%fsc reads Fourier coefficients: transform, correlate, and
        ! return the halves to real space for the evidence setup
        call even%fft()
        call odd%fft()
        call even%fsc(odd, corrs)
        call even%ifft()
        call odd%ifft()
        res = get_resarr(box, smpd)
        call get_resolution(corrs, res, fsc05, fsc0143)
        write(logfhandle,'(A,F8.3,A,F8.3,A)') '>>> POSTPROCESS_NU: HALF-MAP FSC=0.5 ', fsc05, ' A, FSC=0.143 ', fsc0143, ' A'
        if( fsc0143 <= TINY .or. fsc0143 <= 2.*smpd + TINY ) &
            &THROW_HARD('the half-map FSC never falls below 0.143: are the _even_unfil/_odd_unfil halves independent? postprocess_nu')
        if( l_aux )then
            ! the refinement's filter competition: static ladder capped at
            ! fsc/1.5 plus the regularized pair beside the finest rung;
            ! its products (_nu_filt references, _nu_locres), the assignment
            ! table and the matching handoff are exactly what a refinement
            ! iteration would use (2026-09-18)
            write(logfhandle,'(A)') '>>> POSTPROCESS_NU: FILTER COMPETITION WITH THE ML-REGULARIZED PAIR'
            call setup_nu_dmats(even, odd, params%mskdiam, [fsc0143], aux_even, aux_odd, fsc_res=fsc0143)
            call optimize_nu_cutoff_finds()
            call print_nu_filtmap_lowpass_stats()
            raw_lp     = get_nu_filtmap_finest_selected_lp(min_assigned_pct=0.)
            handoff_lp = get_nu_filtmap_finest_selected_lp(min_assigned_pct=0., &
                &min_signal_pct=NU_ALIGN_LP_MIN_SIGNAL_PCT, n_signal=n_signal)
            write(logfhandle,'(A,F6.2,A,F6.2,A)') '>>> NU MATCHING LOW-PASS HANDOFF: ', handoff_lp, &
                &' A (raw finest label ', raw_lp, ' A)'
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
            ! every product of this program carries the PPROC_NU_SUFFIX
            fname = basename(add2fbody(fname_odd, params%ext, PPROC_NU_SUFFIX//'_filt'))
            call vol_odd_nu%write(fname, del_if_exists=.true.)
            fname = basename(add2fbody(fname_even, params%ext, PPROC_NU_SUFFIX//'_filt'))
            call vol_even_nu%write(fname, del_if_exists=.true.)
            call vol_even_nu%add(vol_odd_nu)
            call vol_even_nu%mul(0.5)
            fname = basename(add2fbody(fname_vol, params%ext, PPROC_NU_SUFFIX//'_filt'))
            call vol_even_nu%write(fname, del_if_exists=.true.)
            fname = basename(add2fbody(fname_vol, params%ext, PPROC_NU_SUFFIX//'_locres'))
            call write_nu_local_resolution_map(fname)
            write(logfhandle,'(A)') '>>> POSTPROCESS_NU: WROTE THE '//PPROC_NU_SUFFIX//'_filt REFERENCES AND THE '//&
                &PPROC_NU_SUFFIX//'_locres MAP'
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call aux_odd(1)%kill
            call aux_even(1)%kill
            deallocate(aux_odd, aux_even)
            call cleanup_nu_filter()
            call fname%kill
        endif
        ! frozen evidence from the unregularized half pair, the standard
        ! lifecycle: bank -> compact immutable state (one evidence identity,
        ! no second NU analysis)
        call setup_nu_dmats(even, odd, params%mskdiam, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE, &
            &fsc_res=fsc0143)
        call optimize_nu_cutoff_finds()
        write(logfhandle,'(A,F8.3,A)') '>>> POSTPROCESS_NU: EVIDENCE BANK FINEST MEMBER ', get_nu_filter_bank_finest_lp(), ' A'
        call build_nu_evidence_state(even, odd, evstate)
        call cleanup_nu_filter()
        call assert_nu_evidence_replay_ready(evstate)
        call print_nu_evidence_summary(evstate)
        ! classical shrink-then-sharpen localized by the evidence; the shipped
        ! product is the single sharpened merged volume
        call nu_evidence_sharpen_vol(evstate, even, odd, vol_sharp)
        if( params%outvol .ne. '' )then
            vol_out = params%outvol
        else
            vol_out = basename(add2fbody(fname_vol, params%ext, PPROC_NU_SUFFIX))
        endif
        call vol_sharp%write(vol_out, del_if_exists=.true.)
        call wait_for_closure(vol_out)
        ! destruct
        call even%kill
        call odd%kill
        call vol_sharp%kill
        call vol_out%kill
        call fname_vol%kill
        call fname_even_unfil%kill
        call fname_odd_unfil%kill
        call fname_even%kill
        call fname_odd%kill
        if( allocated(corrs) ) deallocate(corrs)
        if( allocated(res)   ) deallocate(res)
        call simple_end('**** SIMPLE_POSTPROCESS_NU NORMAL STOP ****')
    end subroutine exec_postprocess_nu

end module simple_commanders_postprocess_nu
