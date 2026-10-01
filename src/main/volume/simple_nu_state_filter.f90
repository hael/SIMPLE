!@descr: assembly-owned nonuniform (NU) filtering of one state's half-map pair
! Shared by gridding volassemble and the PCG master: bank from the base pair with the ML pair as optional
! auxiliary member (setup_nu_dmats), the automsk background (density, or valid NU evidence for automsk=nu),
! _nu_filt references (composed from the apply pair when given, pcg_solvent=yes), _nu_locres, and the
! handoff (finest bank member). Contract: doc/policies/NU/nonuniform_filtering_policy.md.
module simple_nu_state_filter
use simple_core_module_api
use simple_image,            only: image
use simple_image_msk,        only: image_msk
use simple_parameters,       only: parameters
use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
    &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
    &NU_DEV_OUTPUT, get_nu_filtmap_finest_selected_lp, NU_ALIGN_LP_MIN_SIGNAL_PCT, &
    &get_nu_filter_bank_finest_lp, &
    &write_nu_local_resolution_map, write_nu_evidence_envmask, &
    &set_nu_evidence_null_shell, set_nu_solvent_envelope
implicit none

public :: nonuniform_filter_state, nu_state_filter_timings, nu_aux_member
private
#include "simple_local_flags.inc"

type nu_state_filter_timings
    real(timer_int_kind) :: envmask = 0.
    real(timer_int_kind) :: filter  = 0.
end type nu_state_filter_timings

contains

    !> The ML-regularized pair is offered as the auxiliary member whenever ml_reg=yes;
    !! setup_nu_dmats decides whether it joins the bank.
    pure logical function nu_aux_member( params ) result( l_use_aux )
        class(parameters), intent(in) :: params
        l_use_aux = params%l_ml_reg
    end function nu_aux_member

    !> NU competition for one state on the base pair; the base and ML pairs are both killed here, the ML pair
    !! offered to the bank only when l_use_aux. res0143 caps the bank and is the auxiliary's resolution;
    !! align_lp = finest bank member. An optional apply pair (kept) receives the label field instead.
    subroutine nonuniform_filter_state( params, state, vol_base_even, vol_base_odd, &
            &vol_aux_even, vol_aux_odd, l_use_aux, res0143, volname, eonames, align_lp, timings, &
            &base_support, l_base_constrained, vol_apply_even, vol_apply_odd )
        class(parameters),            intent(in)    :: params
        integer,                      intent(in)    :: state
        type(image),                  intent(inout) :: vol_base_even, vol_base_odd
        type(image),                  intent(inout) :: vol_aux_even, vol_aux_odd
        logical,                      intent(in)    :: l_use_aux
        real,                         intent(in)    :: res0143
        class(string),                intent(in)    :: volname, eonames(2)
        real,                         intent(out)   :: align_lp
        type(nu_state_filter_timings), optional, intent(inout) :: timings
        class(image),     optional, intent(in)    :: base_support       !< the support that constrained the base pair (PCG)
        logical,          optional, intent(in)    :: l_base_constrained !< the base pair was solved on base_support, not the sphere
        class(image),     optional, intent(in)    :: vol_apply_even, vol_apply_odd !< the pair the label field is applied to
        type(image), allocatable :: nu_aux_even(:), nu_aux_odd(:)
        type(image)              :: vol_even_nu, vol_odd_nu, vol_base_avg, envelope_core, envelope_dilated
        type(image_msk)          :: density_envelope, active_envelope
        type(string)             :: nu_envmask_file
        integer(timer_int_kind)  :: t_filter, t_envmask
        real    :: aux_resolution, bank_cap_res
        logical :: l_constrained, l_nu_envelope_valid, l_apply
        align_lp = 0.
        if( L_BENCH_GLOB ) t_filter = tic()
        l_constrained = .false.
        if( present(l_base_constrained) ) l_constrained = l_base_constrained
        if( l_constrained .and. .not.present(base_support) ) &
            &THROW_HARD('an envelope-constrained base pair must be accompanied by its base support; nonuniform_filter_state')
        if( trim(params%automsk).ne.'no' )then
            ! conservative density envelope of the base pair (automask3D at envmsklp), built before setup
            ! consumes the pair; the core/dilated intermediates are kept only for the Euclidean null shell
            call vol_base_avg%copy(vol_base_even)
            call vol_base_avg%add(vol_base_odd)
            call vol_base_avg%mul(0.5)
            if( l_constrained )then
                call density_envelope%automask3D(params, vol_base_avg, .false., lp_override=params%envmsklp, &
                    &l_report=.false., core=envelope_core, dilated=envelope_dilated)
            else
                call density_envelope%automask3D(params, vol_base_avg, .false., lp_override=params%envmsklp, &
                    &l_report=.false.)
            endif
            call vol_base_avg%kill
            if( trim(params%automsk) == 'nu' ) &
                &call density_envelope%write(string(AUTOMASK_FBODY//int2str_pad(state,2)//MRC_EXT), del_if_exists=.true.)
        endif
        ! candidate bank from the base pair (the static ladder cut at
        ! fsc/NU_BANK_FSC_HEADROOM); the ML pair carries its FSC=0.143
        ! resolution and joins beside the finest rung once that is at or
        ! beyond it (setup_nu_dmats, one rule for every workflow)
        bank_cap_res = res0143
        if( l_use_aux )then
            allocate(nu_aux_even(1), nu_aux_odd(1))
            call nu_aux_even(1)%copy(vol_aux_even)
            call nu_aux_odd(1)%copy(vol_aux_odd)
            aux_resolution = res0143
            call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [aux_resolution], &
                &nu_aux_even, nu_aux_odd, fsc_res=bank_cap_res)
        else
            call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [real ::], &
                &fsc_res=bank_cap_res)
        endif
        if( trim(params%automsk).ne.'no' )then
            ! evidence envelope from the live unaries: background and reference mask under automsk=nu when
            ! valid (density fallback), a diagnostic under automsk=yes; null regime per set_nu_evidence_null_shell
            if( L_BENCH_GLOB ) t_envmask = tic()
            if( l_constrained ) call set_nu_evidence_null_shell(density_envelope, envelope_core, envelope_dilated, base_support)
            nu_envmask_file = string(NU_ENVMASK_FBODY)//int2str_pad(state,2)//string(MRC_EXT)
            l_nu_envelope_valid = .false.
            call write_nu_evidence_envmask(params%nu_msk_sig, params%amsklp, vol_base_even%get_smpd(), &
                &state, nu_envmask_file, mask_out=active_envelope, l_valid=l_nu_envelope_valid)
            if( trim(params%automsk) == 'nu' .and. l_nu_envelope_valid )then
                call set_nu_solvent_envelope(active_envelope, source='nu_evidence_envelope')
                write(logfhandle,'(A,I0)') &
                    &'>>> NU BACKGROUND: COARSEST CANDIDATE OUTSIDE THE NU EVIDENCE ENVELOPE, STATE ', state
            else
                if( trim(params%automsk) == 'nu' )then
                    if( file_exists(nu_envmask_file) ) call del_file(nu_envmask_file)
                    write(logfhandle,'(A,I0,A)') '>>> NU BACKGROUND: STATE ', state, &
                        &', NU evidence envelope unavailable or invalid; using density fallback'
                endif
                call active_envelope%copy(density_envelope)
                call set_nu_solvent_envelope(active_envelope, source='density_envelope')
                write(logfhandle,'(A,I0)') &
                    &'>>> NU BACKGROUND: COARSEST CANDIDATE OUTSIDE THE DENSITY ENVELOPE, STATE ', state
            endif
            call nu_envmask_file%kill
            if( l_constrained )then
                call envelope_core%kill
                call envelope_dilated%kill
            endif
            if( L_BENCH_GLOB .and. present(timings) ) timings%envmask = timings%envmask + toc(t_envmask)
        endif
        ! the auxiliary inputs are copied into the bank; release them
        call cleanup_nu_aux_images()
        call vol_aux_even%kill
        call vol_aux_odd%kill
        call optimize_nu_cutoff_finds()
        call vol_base_even%kill
        call vol_base_odd%kill
        l_apply = present(vol_apply_even) .and. present(vol_apply_odd)
        if( present(vol_apply_even) .neqv. present(vol_apply_odd) ) &
            &THROW_HARD('an apply pair needs both halves; nonuniform_filter_state')
        if( l_apply )then
            call nu_filter_vols(vol_even_nu, vol_odd_nu, vol_apply_even, vol_apply_odd)
            if( params%part == 1 ) write(logfhandle,'(A,I0,A)') '>>> NU REFERENCES: STATE ', state, &
                &', LABEL FIELD OF THE BASE PAIR APPLIED TO THE SOLVENT-PRIOR PAIR'
        else
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
        endif
        if( trim(params%automsk).ne.'no' )then
            ! The _nu_filt matching references carry the active envelope. In
            ! nu mode this is the current evidence mask, with density fallback.
            call active_envelope%apply_3Dmask(vol_even_nu)
            call active_envelope%apply_3Dmask(vol_odd_nu)
            call active_envelope%kill_bimg
            call density_envelope%kill_bimg
            write(logfhandle,'(A,I0,A,A)') '>>> NU REFERENCES: STATE ', state, ', MULTIPLIED BY THE ', &
                &merge('NU EVIDENCE ENVELOPE', 'DENSITY ENVELOPE    ', &
                    &trim(params%automsk) == 'nu' .and. l_nu_envelope_valid)
        endif
        call print_nu_filtmap_lowpass_stats()
        if( NU_DEV_OUTPUT .and. params%part == 1 ) call analyze_filtmap_neighbor_continuity()
        call write_nonuniform_outputs()
        call record_nu_alignment_lowpass_limit()
        call vol_even_nu%kill
        call vol_odd_nu%kill
        call cleanup_nu_filter()
        if( L_BENCH_GLOB .and. present(timings) ) timings%filter = timings%filter + toc(t_filter)

    contains

        subroutine write_nonuniform_outputs()
            type(string) :: eonames_nu(2), volname_nu, locres_name
            eonames_nu(1) = add2fbody(eonames(1), MRC_EXT, NUFILT_SUFFIX)
            eonames_nu(2) = add2fbody(eonames(2), MRC_EXT, NUFILT_SUFFIX)
            volname_nu    = add2fbody(volname,    MRC_EXT, NUFILT_SUFFIX)
            locres_name   = add2fbody(volname,    MRC_EXT, NULOCRES_SUFFIX)
            call vol_even_nu%write(eonames_nu(1), del_if_exists=.true.)
            call vol_odd_nu%write(eonames_nu(2), del_if_exists=.true.)
            call vol_even_nu%add(vol_odd_nu)
            call vol_even_nu%mul(0.5)
            call vol_even_nu%write(volname_nu, del_if_exists=.true.)
            call write_nu_local_resolution_map(locres_name)
            call wait_for_closure(volname_nu)
            call wait_for_closure(locres_name)
            call eonames_nu(1)%kill
            call eonames_nu(2)%kill
            call volname_nu%kill
            call locres_name%kill
        end subroutine write_nonuniform_outputs

        subroutine record_nu_alignment_lowpass_limit()
            real    :: raw_lp, populated_lp
            integer :: n_signal
            ! handoff = finest bank member, independent of which labels won voxels; the population statistics
            ! are diagnostics (doc/policies/NU/nonuniform_filtering_policy.md section 12)
            align_lp     = get_nu_filter_bank_finest_lp()
            raw_lp       = get_nu_filtmap_finest_selected_lp(min_assigned_pct=0.)
            populated_lp = get_nu_filtmap_finest_selected_lp(min_assigned_pct=0., &
                &min_signal_pct=NU_ALIGN_LP_MIN_SIGNAL_PCT, n_signal=n_signal)
            if( params%part == 1 )then
                write(logfhandle,'(A,I0,A,F6.2,A,F6.2,A,F6.2,A)') &
                    &'>>> NU MATCHING LOW-PASS HANDOFF: STATE ', state, ', ', align_lp, &
                    &' A (finest bank member; finest label with 1% of signal voxels ', populated_lp, &
                    &' A, raw finest label ', raw_lp, ' A)'
            endif
        end subroutine record_nu_alignment_lowpass_limit

        subroutine cleanup_nu_aux_images()
            integer :: i
            if( allocated(nu_aux_even) )then
                do i = 1, size(nu_aux_even)
                    call nu_aux_even(i)%kill
                enddo
                deallocate(nu_aux_even)
            endif
            if( allocated(nu_aux_odd) )then
                do i = 1, size(nu_aux_odd)
                    call nu_aux_odd(i)%kill
                enddo
                deallocate(nu_aux_odd)
            endif
        end subroutine cleanup_nu_aux_images

    end subroutine nonuniform_filter_state

end module simple_nu_state_filter
