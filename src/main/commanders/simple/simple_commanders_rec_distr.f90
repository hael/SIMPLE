module simple_commanders_rec_distr
use simple_commanders_api
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody, refine3D_resolution_txt_fbody, &
    &refine3D_state_halfvol_fname, refine3D_state_vol_fname, refine3D_fsc_fname, &
    &refine3D_volassemble_bench_fname, refine3D_trail_rec_fbody, refine3D_trail_rec_fname, &
    &refine3D_trail_rho_fname, refine3D_trail_manifest_fname
implicit none
private
public :: commander_volassemble, filter_pcg_nonuniform_maps
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_volassemble
  contains
    procedure :: execute => exec_volassemble
end type commander_volassemble

type :: restore_timings_t
    real(timer_int_kind) :: reduce_partials                = 0.
    real(timer_int_kind) :: sum_eos                        = 0.
    real(timer_int_kind) :: restore_eos_and_write_fsc      = 0.
    real(timer_int_kind) :: restore_merged_volume          = 0.
    real(timer_int_kind) :: trail_blend_accums             = 0.
    real(timer_int_kind) :: trail_restored_halves          = 0.
end type restore_timings_t

contains

    !> Reduce one state's Cartesian partial reconstructions and restore dense
    !> even, odd, and merged volumes. On return build%vol/build%vol2 contain the
    !> restored half-volumes needed by automask, while vol_nu_base_even/odd and
    !> optional vol_nu_aux_even/odd contain the static-bank nonuniform-filter
    !> auxiliary inputs before even/odd low-resolution insertion.
    subroutine restore_state_from_parts( params, build, cline, even_rec, odd_rec, read_even_rec, read_odd_rec, &
        &sum_rec, state, numlen_part, &
        &update_frac_trail_rec, realized_update_frac, vol_prev_even, vol_prev_odd, vol_merged, &
        &vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd, &
        &volname, eonames, res05, res0143, timings )
        use simple_reconstructor, only: reconstructor
        use simple_reconstructor_eo, only: gridding_pair_diagnostics, calculate_gridding_pair_fsc, &
            &read_gridding_pair_accumulators, restore_gridding_pair, write_gridding_pair_accumulators, &
            &write_gridding_pair_diagnostics
        type(parameters),       intent(in)    :: params
        type(builder),          intent(inout) :: build
        class(cmdline),         intent(in)    :: cline
        type(reconstructor),    intent(inout) :: even_rec, odd_rec, read_even_rec, read_odd_rec, sum_rec
        integer,                intent(in)    :: state, numlen_part
        real,                   intent(in)    :: update_frac_trail_rec !< applied map-update weight u (ufrac_trec override or realized)
        real,                   intent(in)    :: realized_update_frac  !< realized fraction f that produced the current partials
        type(image),            intent(inout) :: vol_prev_even, vol_prev_odd, vol_merged
        type(image),            intent(inout) :: vol_nu_base_even, vol_nu_base_odd
        type(image),            intent(inout) :: vol_nu_aux_even, vol_nu_aux_odd
        type(string),           intent(inout) :: volname, eonames(2)
        real,                   intent(out)   :: res05, res0143
        type(restore_timings_t), intent(inout) :: timings
        type(string) :: volname_prev, volname_prev_even, volname_prev_odd
        type(string) :: fsc_txt_file, trail_fbody
        type(gridding_pair_diagnostics) :: pair_diagnostics
        real, allocatable :: fsc(:)
        real    :: weight_prev, bootstrap_cfar
        integer :: ldim(3), trail_chain_gen
        logical :: l_trail_chain
        integer(timer_int_kind) :: t_reduce_partials, t_restore_eos, t_restore_merged, t_sum_eos, t_trail
        integer(timer_int_kind) :: t_trail_blend
        call reduce_partials()
        call set_state_filenames()
        ! Trailing happens in the accumulator domain: the persistent chain of
        ! e/o Fourier sums + sampling densities is decayed by (1 - update
        ! fraction) and the current partial sums are added, so every Fourier
        ! component is weighted by its accumulated sampling density. All
        ! downstream restoration (FSC, halfmaps, merged volume, nonuniform
        ! filter inputs) then consumes the blended statistics. Only when the
        ! chain does not exist yet (bootstrap) do we fall back to the legacy
        ! previous-halfmap volume blend below, while seeding the chain.
        call blend_trailing_accumulators()
        call sum_eos_before_density_correction_if_needed()
        call restore_eos_and_write_fsc()
        call sum_eos_after_density_correction_if_needed()
        call capture_nonuniform_source_halves()
        call restore_merged_volume()
        ! Keep restored half-volumes current in build%vol/build%vol2 for
        ! automasking and optional trailing. Low-resolution even/odd blending is
        ! a registration-reference trick and is applied only during
        ! reprojection-model preparation, never to these on-disk halfmaps.
        call build%vol%read(eonames(1))
        call build%vol2%read(eonames(2))
        call trail_restored_halves_if_needed()
        call cleanup_restore_state()

    contains

        subroutine reduce_partials()
            integer :: part
            if( L_BENCH_GLOB ) t_reduce_partials = tic()
            call even_rec%reset
            call odd_rec%reset
            do part = 1, params%nparts
                call read_gridding_pair_accumulators(params, read_even_rec, read_odd_rec, &
                    &refine3D_partial_rec_fbody(state, part, numlen_part))
                call even_rec%sum_reduce(read_even_rec)
                call odd_rec%sum_reduce(read_odd_rec)
            enddo
            if( L_BENCH_GLOB ) timings%reduce_partials = timings%reduce_partials + toc(t_reduce_partials)
        end subroutine reduce_partials

        subroutine set_state_filenames()
            volname    = refine3D_state_vol_fname(state)
            eonames(1) = refine3D_state_halfvol_fname(state, 'even')
            eonames(2) = refine3D_state_halfvol_fname(state, 'odd')
        end subroutine set_state_filenames

        subroutine blend_trailing_accumulators()
            real :: cur_scale
            l_trail_chain   = .false.
            trail_chain_gen = 0
            trail_fbody     = refine3D_trail_rec_fbody(state)
            if( .not. params%l_trail_rec )then
                ! Full-reconstruction producer contract: a stage-boundary
                ! reconstruct3D can seed the chain with full-dataset weight via
                ! the internal trail_seed handshake, so the consuming trailing
                ! stage starts from complete statistics instead of warming up.
                if( cline%defined('trail_seed') )then
                    if( cline%get_carg('trail_seed') .eq. 'yes' )then
                        call write_trail_chain_set()
                        write(logfhandle,'(A,I0)') &
                            &'>>> VOLASSEMBLE: WROTE FULL-WEIGHT TRAILING CHAIN SEED, STATE ', state
                    endif
                endif
                return
            endif
            l_trail_chain = validate_trail_chain()
            if( realized_update_frac < 0.001 )then
                ! nothing meaningful was sampled for this state; leave the chain
                ! untouched and let the legacy path govern this iteration
                THROW_WARN('realized update fraction ~0; trailing chain left untouched')
                l_trail_chain = .false.
                return
            endif
            weight_prev = 1.0 - update_frac_trail_rec
            if( l_trail_chain .and. weight_prev > 0.01 )then
                if( L_BENCH_GLOB ) t_trail_blend = tic()
                ! ufrac_trec contract: the applied fraction u must be the restored
                ! current-map coefficient even when it differs from the realized
                ! fraction f that produced the partials. Scaling the current
                ! accumulators by u/f achieves that and keeps the chain at full
                ! sampling mass: (u/f)*(f*D) + (1-u)*D = D.
                cur_scale = update_frac_trail_rec / realized_update_frac
                if( abs(cur_scale - 1.0) > 0.001 )then
                    call even_rec%apply_weight_sums(cur_scale)
                    call odd_rec%apply_weight_sums(cur_scale)
                endif
                call read_gridding_pair_accumulators(params, read_even_rec, read_odd_rec, trail_fbody, &
                    &required=.true.)
                call read_even_rec%apply_weight_sums(weight_prev)
                call read_odd_rec%apply_weight_sums(weight_prev)
                call even_rec%sum_reduce(read_even_rec)
                call odd_rec%sum_reduce(read_odd_rec)
                call write_trail_chain_set()
                write(logfhandle,'(A,I0,A,F8.4,A,F8.4)') &
                    &'>>> VOLASSEMBLE: TRAILING ACCUMULATOR BLEND, STATE ', state, &
                    &', PREVIOUS-CHAIN WEIGHT ', weight_prev, ', CURRENT SCALE ', cur_scale
                if( L_BENCH_GLOB ) timings%trail_blend_accums = &
                    timings%trail_blend_accums + toc(t_trail_blend)
            else
                ! No blend this iteration: either the chain does not exist yet
                ! (bootstrap; the legacy previous-halfmap volume blend produces
                ! this iteration's outputs) or the applied fraction is ~1 (full
                ! replacement of the model). Either way the persisted chain must
                ! represent full-dataset sampling mass: fractional partials are
                ! scaled by 1/f for the write and restored afterwards. Without
                ! this, a chain seeded at fractional mass f makes the next
                ! iteration's effective update weight f/(f + (1-f)*f), far above
                ! the requested fraction.
                call even_rec%apply_weight_sums(1.0 / realized_update_frac)
                call odd_rec%apply_weight_sums(1.0 / realized_update_frac)
                call write_trail_chain_set()
                call even_rec%apply_weight_sums(realized_update_frac)
                call odd_rec%apply_weight_sums(realized_update_frac)
                if( .not. l_trail_chain )then
                    write(logfhandle,'(A,I0,A)') &
                        &'>>> VOLASSEMBLE: SEEDED FULL-MASS TRAILING CHAIN, STATE ', state, &
                        &'; USING LEGACY PREVIOUS-HALFMAP BLEND THIS ITERATION'
                endif
            endif
        end subroutine blend_trailing_accumulators

        function trail_chain_component( ifile ) result( fname )
            integer, intent(in) :: ifile
            type(string) :: fname
            select case(ifile)
                case(1)
                    fname = refine3D_trail_rec_fname(state, 'even')
                case(2)
                    fname = refine3D_trail_rho_fname(state, 'even')
                case(3)
                    fname = refine3D_trail_rec_fname(state, 'odd')
                case(4)
                    fname = refine3D_trail_rho_fname(state, 'odd')
                case DEFAULT
                    THROW_HARD('invalid trailing chain component index')
            end select
        end function trail_chain_component

        integer(kind=8) function trail_chain_file_size( fname ) result( sz )
            class(string), intent(in) :: fname
            sz = -1
            if( file_exists(fname) ) inquire(file=fname%to_char(), size=sz)
        end function trail_chain_file_size

        !> Remove the complete artifact set (manifest + four accumulator files)
        !! so an invalid chain can never contribute components to a later one.
        subroutine discard_trail_chain_set( reason )
            character(len=*), intent(in), optional :: reason
            type(string) :: fname
            logical      :: l_any
            integer      :: ifile
            l_any = .false.
            do ifile = 1, 4
                fname = trail_chain_component(ifile)
                if( file_exists(fname) )then
                    l_any = .true.
                    call del_file(fname)
                endif
            enddo
            fname = refine3D_trail_manifest_fname(state)
            if( file_exists(fname) )then
                l_any = .true.
                call del_file(fname)
            endif
            call fname%kill
            if( l_any .and. present(reason) )then
                write(logfhandle,'(A,I0,A)') &
                    &'>>> VOLASSEMBLE: DISCARDING TRAILING CHAIN, STATE ', state, ': '//trim(reason)
            endif
        end subroutine discard_trail_chain_set

        !> Persist the chain as one artifact set. The manifest is deleted first
        !! and rewritten last with the byte size of every component, so an
        !! interrupted write leaves no valid manifest and the mixed-generation
        !! remnants are rejected and discarded by the next validation.
        subroutine write_trail_chain_set()
            type(string)    :: manifest
            integer(kind=8) :: sizes(4)
            integer         :: funit, io_stat, ifile
            manifest = refine3D_trail_manifest_fname(state)
            call del_file(manifest)
            call write_gridding_pair_accumulators(even_rec, odd_rec, trail_fbody)
            do ifile = 1, 4
                sizes(ifile) = trail_chain_file_size(trail_chain_component(ifile))
            enddo
            call fopen(funit, file=manifest, status='REPLACE', action='WRITE', iostat=io_stat)
            if( io_stat == 0 )then
                write(funit,*) params%box_crop, params%smpd_crop, build%spproj_field%get_noris(), &
                    &params%nstates, state, trail_chain_gen + 1, sizes
                call fclose(funit)
            else
                THROW_WARN('failed to write trailing chain manifest; chain will be re-seeded next iteration')
            endif
            call manifest%kill
        end subroutine write_trail_chain_set

        !> The chain is one validated artifact set: the manifest must exist,
        !! parse, match the current project population and state layout, match
        !! every component's byte size, and describe a grid that is not larger
        !! than the current one with the same physical extent. Smaller previous
        !! grids are zero-padded by the reader (autoscale ramp). Any failure
        !! discards the complete set and re-seeds instead of failing the run.
        logical function validate_trail_chain() result( l_valid )
            type(string)    :: manifest
            integer(kind=8) :: sizes_chain(4)
            integer         :: funit, io_stat, ifile
            integer         :: box_chain, nptcls_chain, nstates_chain, state_chain, gen_chain
            real            :: smpd_chain, extent_chain, extent_cur
            l_valid  = .false.
            manifest = refine3D_trail_manifest_fname(state)
            if( .not. file_exists(manifest) )then
                call discard_trail_chain_set() ! remove orphaned components quietly
                call manifest%kill
                return
            endif
            call fopen(funit, file=manifest, status='OLD', action='READ', iostat=io_stat)
            if( io_stat /= 0 )then
                call discard_trail_chain_set('unreadable manifest')
                call manifest%kill
                return
            endif
            read(funit,*,iostat=io_stat) box_chain, smpd_chain, nptcls_chain, nstates_chain, &
                &state_chain, gen_chain, sizes_chain
            call fclose(funit)
            call manifest%kill
            if( io_stat /= 0 )then
                call discard_trail_chain_set('corrupt manifest')
                return
            endif
            trail_chain_gen = gen_chain
            if( nptcls_chain /= build%spproj_field%get_noris() )then
                call discard_trail_chain_set('particle population mismatch')
                return
            endif
            if( nstates_chain /= params%nstates .or. state_chain /= state )then
                call discard_trail_chain_set('state layout mismatch')
                return
            endif
            if( box_chain > params%box_crop )then
                call discard_trail_chain_set('larger grid than current reconstruction')
                return
            endif
            extent_chain = real(box_chain)      * smpd_chain
            extent_cur   = real(params%box_crop) * params%smpd_crop
            if( abs(extent_chain - extent_cur) > 0.01 * extent_cur )then
                call discard_trail_chain_set('physical extent mismatch')
                return
            endif
            do ifile = 1, 4
                if( trail_chain_file_size(trail_chain_component(ifile)) /= sizes_chain(ifile) )then
                    call discard_trail_chain_set('component size mismatch (mixed generations)')
                    return
                endif
            enddo
            l_valid = .true.
        end function validate_trail_chain

        subroutine sum_eos_before_density_correction_if_needed()
            if( params%l_ml_reg ) return
            call sum_pair_into_sum_rec()
        end subroutine sum_eos_before_density_correction_if_needed

        !> single owner of the merged-accumulator summation so the ml_reg and
        !! non-ml_reg paths cannot diverge
        subroutine sum_pair_into_sum_rec()
            if( L_BENCH_GLOB ) t_sum_eos = tic()
            call sum_rec%reset
            call sum_rec%sum_reduce(even_rec)
            call sum_rec%sum_reduce(odd_rec)
            if( L_BENCH_GLOB ) timings%sum_eos = timings%sum_eos + toc(t_sum_eos)
        end subroutine sum_pair_into_sum_rec

        subroutine restore_eos_and_write_fsc()
            use simple_fsc, only: fsc_area_score_result
            type(fsc_area_score_result) :: cones_fsc
            if( L_BENCH_GLOB ) t_restore_eos = tic()
            if( params%l_trail_rec .and. .not. l_trail_chain )then
                ! Bootstrap iteration: the chain has no information yet, so the
                ! previous halfmaps provide the FSC prior for regularization,
                ! exactly as the legacy volume-domain trailing did.
                call read_previous_halfmaps()
                if( allocated(fsc) ) deallocate(fsc)
                if( params%conical_fsc == 'yes' )then
                    call cones_fsc%new(vol_prev_even, 256, 20., 0.143, 1)
                    call calculate_gridding_pair_fsc(params, vol_prev_even, vol_prev_odd, &
                        &state, fsc, bootstrap_cfar, cones=cones_fsc)
                    call restore_gridding_pair(params, even_rec, odd_rec, state, &
                        &eonames(1), eonames(2), pair_diagnostics, fsc_in=fsc, &
                        &cfar_in=bootstrap_cfar, cones_in=cones_fsc)
                    call cones_fsc%kill
                else
                    call calculate_gridding_pair_fsc(params, vol_prev_even, vol_prev_odd, &
                        &state, fsc, bootstrap_cfar)
                    call restore_gridding_pair(params, even_rec, odd_rec, state, &
                        &eonames(1), eonames(2), pair_diagnostics, fsc_in=fsc, cfar_in=bootstrap_cfar)
                endif
            else
                ! With an accumulator chain the blended sums already contain the
                ! trailed statistics, so the FSC is estimated post-blend from
                ! the restored halves and describes the artifact written to disk.
                call restore_gridding_pair(params, even_rec, odd_rec, state, &
                    &eonames(1), eonames(2), pair_diagnostics)
            endif
            res05      = pair_diagnostics%res_fsc05
            res0143    = pair_diagnostics%res_fsc0143
            fsc_txt_file = resolve_fsc_txt_fname()
            call write_gridding_pair_diagnostics(pair_diagnostics, params%box_crop, params%smpd_crop, fsc_txt_file)
            if( L_BENCH_GLOB ) timings%restore_eos_and_write_fsc = &
                timings%restore_eos_and_write_fsc + toc(t_restore_eos)
        end subroutine restore_eos_and_write_fsc

        subroutine read_previous_halfmaps()
            if( .not. cline%defined('vol'//int2str(state)) )then
                THROW_HARD('vol'//int2str(state)//' required in volassemble cline when trail_rec==yes')
            endif
            volname_prev      = cline%get_carg('vol'//int2str(state))
            volname_prev_even = add2fbody(volname_prev, MRC_EXT, '_even')
            volname_prev_odd  = add2fbody(volname_prev, MRC_EXT, '_odd')
            if( .not. file_exists(volname_prev_even) )then
                THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
            endif
            if( .not. file_exists(volname_prev_odd) )then
                THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
            endif
            call vol_prev_even%read_and_crop(volname_prev_even, params%smpd, params%box_crop, params%smpd_crop)
            call vol_prev_odd%read_and_crop( volname_prev_odd,  params%smpd, params%box_crop, params%smpd_crop)
        end subroutine read_previous_halfmaps

        function resolve_fsc_txt_fname() result( fname )
            type(string) :: fname
            if( cline%defined('outfile') )then
                fname = resolution_outfile_fbody()
            else if( cline%defined('which_iter') )then
                fname = refine3D_resolution_txt_fbody(state, params%which_iter)
            else
                fname = refine3D_resolution_txt_fbody(state)
            endif
        end function resolve_fsc_txt_fname

        function resolution_outfile_fbody() result( fname )
            type(string) :: fname, ext
            fname = params%outfile
            ext   = fname2ext(fname)
            select case(ext%to_char())
                case('txt','simple')
                    fname = get_fbody(fname, ext)
            end select
            fname = fname//'_STATE'//int2str_pad(state,2)
            call ext%kill
        end function resolution_outfile_fbody

        subroutine sum_eos_after_density_correction_if_needed()
            if( .not. params%l_ml_reg ) return
            call sum_pair_into_sum_rec()
        end subroutine sum_eos_after_density_correction_if_needed

        subroutine capture_nonuniform_source_halves()
            if( .not. params%l_nonuniform ) return
            ldim = [params%box_crop, params%box_crop, params%box_crop]
            call vol_nu_base_even%new(ldim, params%smpd_crop)
            call vol_nu_base_odd%new( ldim, params%smpd_crop)
            if( params%l_ml_reg )then
                call vol_nu_base_even%read(add2fbody(eonames(1), MRC_EXT, '_unfil'))
                call vol_nu_base_odd%read( add2fbody(eonames(2), MRC_EXT, '_unfil'))
                if( use_static_nu_aux_replacement() )then
                    call vol_nu_aux_even%new(ldim, params%smpd_crop)
                    call vol_nu_aux_odd%new( ldim, params%smpd_crop)
                    call vol_nu_aux_even%read(eonames(1))
                    call vol_nu_aux_odd%read( eonames(2))
                endif
            else
                call vol_nu_base_even%read(eonames(1))
                call vol_nu_base_odd%read( eonames(2))
            endif
            ! the half-maps on disk are deapodized by the one-half gridding backend; nothing to apply here
        end subroutine capture_nonuniform_source_halves

        subroutine restore_merged_volume()
            if( L_BENCH_GLOB ) t_restore_merged = tic()
            if( L_VERBOSE_GLOB )then
                write(logfhandle,'(A)') '>>> SAMPLING DENSITY (RHO) CORRECTION & WIENER NORMALIZATION'
            endif
            call sum_rec%restore_final(build%vol)
            call build%vol%fft
            call build%vol%ifft
            call build%vol%write(volname, del_if_exists=.true.)
            call wait_for_closure(volname)
            if( L_BENCH_GLOB ) timings%restore_merged_volume = &
                timings%restore_merged_volume + toc(t_restore_merged)
        end subroutine restore_merged_volume

        !> Legacy volume-domain blend, retained ONLY for the bootstrap iteration
        !! that seeds the accumulator chain: it keeps the output halves anchored
        !! to the previous model exactly as before. Chain iterations never enter
        !! here — their halves, merged volume and NU inputs are already trailed.
        subroutine trail_restored_halves_if_needed()
            if( .not. params%l_trail_rec ) return
            if( l_trail_chain ) return
            if( update_frac_trail_rec >= 0.99 ) return
            if( L_BENCH_GLOB ) t_trail = tic()
            weight_prev = 1. - update_frac_trail_rec
            call vol_prev_even%mul(weight_prev)
            call vol_prev_odd%mul(weight_prev)
            call build%vol%mul(update_frac_trail_rec)
            call build%vol2%mul(update_frac_trail_rec)
            call build%vol%add(vol_prev_even)
            call build%vol2%add(vol_prev_odd)
            if( params%l_nonuniform )then
                call vol_nu_base_even%mul(update_frac_trail_rec)
                call vol_nu_base_odd%mul(update_frac_trail_rec)
                call vol_nu_base_even%add(vol_prev_even)
                call vol_nu_base_odd%add(vol_prev_odd)
                if( use_static_nu_aux_replacement() )then
                    call vol_nu_aux_even%mul(update_frac_trail_rec)
                    call vol_nu_aux_odd%mul(update_frac_trail_rec)
                    call vol_nu_aux_even%add(vol_prev_even)
                    call vol_nu_aux_odd%add(vol_prev_odd)
                endif
            endif
            call build%vol%write(eonames(1))
            call build%vol2%write(eonames(2))
            if( params%l_lpset )then
                call vol_merged%copy(build%vol)
                call vol_merged%add(build%vol2)
                call vol_merged%mul(0.5)
                call vol_merged%write(volname, del_if_exists=.true.)
                call wait_for_closure(volname)
            endif
            if( L_BENCH_GLOB ) timings%trail_restored_halves = timings%trail_restored_halves + toc(t_trail)
        end subroutine trail_restored_halves_if_needed

        logical function use_static_nu_aux_replacement() result(l_use_aux)
            l_use_aux = params%l_ml_reg .and. .not. params%l_nu_refine
        end function use_static_nu_aux_replacement

        subroutine cleanup_restore_state()
            call vol_prev_even%kill
            call vol_prev_odd%kill
            call fsc_txt_file%kill
            call trail_fbody%kill
            call volname_prev%kill
            call volname_prev_even%kill
            call volname_prev_odd%kill
            call pair_diagnostics%kill
            if( allocated(fsc) ) deallocate(fsc)
        end subroutine cleanup_restore_state

    end subroutine restore_state_from_parts

    !> Apply the volume-assembly nonuniform-filtering contract to half-maps
    !! restored directly by the PCG backend. Raw (B,D) trailing, including its
    !! chain/bootstrap policy, remains wholly owned by the PCG master; this
    !! consumer must never blend solved maps. PCG maps are already corrected in
    !! the solved image domain, so this path deliberately omits the gridding
    !! instrument-function correction performed by commander_volassemble.
    subroutine filter_pcg_nonuniform_maps( params, build, l_trail_bootstrap, nu_replay_lps )
        use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
            &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
            &set_nu_filter_report, NU_DEV_OUTPUT, get_nu_filtmap_finest_selected_lp, &
            &write_nu_local_resolution_map, nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, &
            &print_nu_envmask_stats, NU_ENVMASK_BETA, NU_ENVMASK_DENS_WEIGHT, NU_ENVMASK_RELATIVE, &
            &NU_ENVMASK_MINVOL_FRAC, NU_ENVMASK_GROW_A, NU_ENVMASK_EDGE_A
        use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess, &
            &NU_ENVMASK_ACTION_REGENERATE
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        logical,          intent(in)    :: l_trail_bootstrap(:)
        real, optional,   intent(in)    :: nu_replay_lps(:)
        type(image)                   :: vol_base_even, vol_base_odd, vol_aux_even, vol_aux_odd
        type(image)                   :: vol_even_nu, vol_odd_nu
        type(image), allocatable      :: nu_aux_even(:), nu_aux_odd(:)
        type(image_msk)               :: nu_envmask
        type(vol_pproc_plan)          :: pp_plan
        type(nu_envmask_params)       :: envp
        type(nu_envmask_stats)        :: envstats
        type(string)                  :: eonames(2), volname, eonames_nu(2), volname_nu, locres_name
        type(string)                  :: fsc_fname
        real, allocatable             :: fsc(:), res(:), nu_align_lps(:)
        integer, allocatable          :: state_pops(:)
        logical, allocatable          :: l_env(:,:,:)
        logical, allocatable          :: l_included(:)
        integer                       :: state, i, n_ccs, n_ccs_kept, grow_px, edge_px, ldim(3)
        real                          :: fsc05, fsc0143, aux_resolution, align_lp, selected_lp

        if( .not. params%l_nonuniform ) return
        if( params%l_nu_refine ) THROW_HARD('PCG nonuniform filtering does not yet support nu_refine=yes')
        if( size(l_trail_bootstrap) /= params%nstates ) &
            &THROW_HARD('PCG nonuniform trailing-bootstrap state input has invalid size')
        ! Policy (2026-08-28, pcg_priors.md): with the Q_NU replay active the
        ! shipped maps are already locally regularized in-solve, so the
        ! post-hoc NU filter is redundant -- no second NU analysis, no
        ! _nu_filt/_nu_locres products, no evidence envelope. Only the LP-set
        ! matching handoff survives, derived from the frozen replay evidence
        ! (the finest evidenced local cutoff per state) instead of a second
        ! filter-bank optimization.
        if( params%l_ml_reg .and. params%pcg_nu_lambda_rel > 0.0 )then
            if( .not. present(nu_replay_lps) ) &
                &THROW_HARD('NU replay is active but no evidence-derived LP handoff was supplied')
            if( size(nu_replay_lps) /= params%nstates ) &
                &THROW_HARD('NU replay LP handoff has invalid size')
            allocate(state_pops(params%nstates), l_included(params%nstates))
            do state = 1, params%nstates
                state_pops(state) = build%spproj_field%get_pop(state, 'state')
            enddo
            l_included = (state_pops > 0) .and. (nu_replay_lps > TINY)
            if( any(l_included) )then
                align_lp = minval(nu_replay_lps, mask=l_included)
                write(logfhandle,'(A,F8.3,A)') &
                    &'>>> PCG NU REPLAY: post-hoc NU filtering skipped (Q_NU in-solve); '//&
                    &'evidence-derived matching low-pass ', align_lp, ' A'
                call build%spproj_field%set_all2single('lp', align_lp)
                call build%spproj%write_segment_inside(params%oritype, params%projfile)
            else
                write(logfhandle,'(A)') &
                    &'>>> PCG NU REPLAY: post-hoc NU filtering skipped (Q_NU in-solve); no evidenced cutoff to hand off'
            endif
            deallocate(state_pops, l_included)
            return
        endif
        call set_nu_filter_report(params%part == 1)
        ldim = [params%box_crop, params%box_crop, params%box_crop]
        res  = get_resarr(params%box_crop, params%smpd_crop)
        allocate(state_pops(params%nstates), nu_align_lps(params%nstates), l_included(params%nstates))
        state_pops  = 0
        nu_align_lps = 0.
        do state = 1, params%nstates
            state_pops(state) = build%spproj_field%get_pop(state, 'state')
            if( state_pops(state) < 1 ) cycle

            eonames(1) = refine3D_state_halfvol_fname(state, 'even')
            eonames(2) = refine3D_state_halfvol_fname(state, 'odd')
            volname    = refine3D_state_vol_fname(state)
            call vol_base_even%new(ldim, params%smpd_crop)
            call vol_base_odd%new( ldim, params%smpd_crop)
            if( params%l_ml_reg .and. .not. l_trail_bootstrap(state) )then
                call vol_base_even%read(refine3D_state_halfvol_fname(state, 'even', unfil=.true.))
                call vol_base_odd%read( refine3D_state_halfvol_fname(state, 'odd',  unfil=.true.))
                call vol_aux_even%new(ldim, params%smpd_crop)
                call vol_aux_odd%new( ldim, params%smpd_crop)
                call vol_aux_even%read(eonames(1))
                call vol_aux_odd%read( eonames(2))
                fsc_fname = refine3D_fsc_fname(state)
                fsc       = file2rarr(fsc_fname)
                call get_resolution(fsc, res, fsc05, fsc0143)
                aux_resolution = fsc0143
                if( params%l_lpset .and. params%lp > TINY ) aux_resolution = min(aux_resolution, params%lp)
                allocate(nu_aux_even(1), nu_aux_odd(1))
                call nu_aux_even(1)%copy(vol_aux_even)
                call nu_aux_odd(1)%copy( vol_aux_odd)
                call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [aux_resolution], &
                    &nu_aux_even, nu_aux_odd)
            else
                if( params%l_ml_reg .and. l_trail_bootstrap(state) )then
                    write(logfhandle,'(A,I0,A)') '>>> PCG NU: STATE ', state, &
                        &' IS A CHAINLESS TRAILING BOOTSTRAP; FILTERING THE COHERENT STANDARD PAIR WITHOUT AN ML AUXILIARY'
                endif
                call vol_base_even%read(eonames(1))
                call vol_base_odd%read( eonames(2))
                call setup_nu_dmats(vol_base_even, vol_base_odd, params%mskdiam, [real ::])
            endif
            call cleanup_aux_images()
            call vol_aux_even%kill
            call vol_aux_odd%kill
            call optimize_nu_cutoff_finds()

            call plan_state_postprocess(params, state, params%which_iter, pp_plan)
            if( pp_plan%l_nu_envmask_incompatible )then
                write(logfhandle,'(A,1X,A)') &
                    &'>>> Existing NU evidence envelope incompatible with current box/sampling, regenerating:', &
                    &pp_plan%nu_envmask_file%to_char()
            endif
            if( pp_plan%nu_envmask_action == NU_ENVMASK_ACTION_REGENERATE )then
                envp%nsigma      = params%nu_msk_sig
                envp%beta        = NU_ENVMASK_BETA
                envp%dens_weight = NU_ENVMASK_DENS_WEIGHT
                envp%lp_smooth   = params%amsklp
                envp%l_relative  = NU_ENVMASK_RELATIVE
                call nu_evidence_envelope(envp, l_env, envstats)
                call print_nu_envmask_stats(envstats)
                grow_px = max(1, nint(NU_ENVMASK_GROW_A / params%smpd_crop))
                edge_px = max(1, nint(NU_ENVMASK_EDGE_A / params%smpd_crop))
                call nu_envmask%kill_bimg
                call nu_envmask%envmask3D_from_lmask(l_env, params%smpd_crop, grow_px, edge_px, &
                    &NU_ENVMASK_MINVOL_FRAC, .true., n_ccs, n_ccs_kept)
                ! an empty evidence field returns without constructing the mask
                ! image; writing it would abort on invalid MRC dimensions. Skip
                ! the write -- every envelope consumer handles absence.
                if( n_ccs_kept < 1 )then
                    THROW_WARN('NU evidence envelope is empty; no envelope mask written this iteration')
                else
                    call nu_envmask%write(pp_plan%nu_envmask_file, del_if_exists=.true.)
                    call wait_for_closure(pp_plan%nu_envmask_file)
                endif
                write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') &
                    &'>>> NU ENVELOPE OCCUPANCY: STATE ', state, ', SUPPORT FRACTION ', &
                    &envstats%pct_signal, ' %, COMPONENTS KEPT ', n_ccs_kept, ' OF ', n_ccs
                if( allocated(l_env) ) deallocate(l_env)
            endif

            call vol_base_even%kill
            call vol_base_odd%kill
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
            call print_nu_filtmap_lowpass_stats()
            if( NU_DEV_OUTPUT .and. params%part == 1 ) call analyze_filtmap_neighbor_continuity()
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
            selected_lp = get_nu_filtmap_finest_selected_lp()
            if( selected_lp > TINY ) nu_align_lps(state) = selected_lp

            call vol_even_nu%kill
            call vol_odd_nu%kill
            call cleanup_nu_filter()
            call nu_envmask%kill_bimg
            call pp_plan%nu_envmask_file%kill
            call eonames(1)%kill
            call eonames(2)%kill
            call eonames_nu(1)%kill
            call eonames_nu(2)%kill
            call volname%kill
            call volname_nu%kill
            call locres_name%kill
            call fsc_fname%kill
            if( allocated(fsc) ) deallocate(fsc)
        enddo

        l_included = (state_pops > 0) .and. (nu_align_lps > TINY)
        if( any(l_included) )then
            align_lp = minval(nu_align_lps, mask=l_included)
            call build%spproj_field%set_all2single('lp', align_lp)
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        call cleanup_aux_images()
        call vol_base_even%kill
        call vol_base_odd%kill
        call vol_aux_even%kill
        call vol_aux_odd%kill
        call vol_even_nu%kill
        call vol_odd_nu%kill
        call cleanup_nu_filter()
        call nu_envmask%kill_bimg
        call pp_plan%nu_envmask_file%kill
        if( allocated(fsc) ) deallocate(fsc)
        if( allocated(res) ) deallocate(res)
        deallocate(state_pops, nu_align_lps, l_included)

    contains

        subroutine cleanup_aux_images()
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
        end subroutine cleanup_aux_images

    end subroutine filter_pcg_nonuniform_maps

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor,    only: reconstructor
        use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
            &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
            &set_nu_filter_report, NU_DEV_OUTPUT, &
            &extend_nu_filter_highres_shell_next, refine_nu_extension_filtmap_ordered_labels, &
            &nu_highres_extension_stats, get_nu_filtmap_finest_selected_lp, &
            &get_nu_filtmap_highres_shell_depth, write_nu_local_resolution_map, &
            &nu_envmask_params, nu_envmask_stats, nu_evidence_envelope, print_nu_envmask_stats, &
            &NU_ENVMASK_BETA, NU_ENVMASK_DENS_WEIGHT, NU_ENVMASK_RELATIVE, NU_ENVMASK_MINVOL_FRAC, &
            &NU_ENVMASK_GROW_A, NU_ENVMASK_EDGE_A
        use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess, &
            &NU_ENVMASK_ACTION_REGENERATE
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters), target      :: params
        type(builder)                 :: build
        type(reconstructor)           :: even_rec, odd_rec, read_even_rec, read_odd_rec, sum_rec
        type(image)                   :: vol_prev_even, vol_prev_odd, vol_merged
        type(image)                   :: vol_even_nu, vol_odd_nu
        type(image)                   :: vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd
        type(image), allocatable      :: nu_aux_even(:), nu_aux_odd(:)
        type(image_msk)               :: nu_envmask
        type(string)                  :: volname, eonames(2)
        type(restore_timings_t)       :: restore_timings
        type(vol_pproc_plan)          :: pp_plan
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: state_pops(:)
        logical, allocatable          :: l_state_dropped(:)
        real, allocatable             :: res0143s(:)
        real, allocatable             :: nu_align_lps(:)
        real, allocatable             :: update_frac_trail_recs(:), realized_update_fracs(:)
        real                          :: res05
        integer                       :: state, numlen_part
        integer(timer_int_kind)       :: t_nu_envmask, t_nonuniform_filter, t_tot
        integer(timer_int_kind)       :: t_init_context, t_trail_frac, t_upd_proj, t_cleanup
        real(timer_int_kind)          :: rt_reduce_partials, rt_sum_eos
        real(timer_int_kind)          :: rt_restore_eos_and_write_fsc, rt_restore_merged_volume
        real(timer_int_kind)          :: rt_trail_blend_accums, rt_trail_restored_halves
        real(timer_int_kind)          :: rt_nu_envmask, rt_nonuniform_filter, rt_tot
        real(timer_int_kind)          :: rt_init_context, rt_trail_frac, rt_gridcorr, rt_upd_proj, rt_cleanup
        call initialize_bench_timers()
        if( L_BENCH_GLOB ) t_init_context = tic()
        call initialize_context()
        if( L_BENCH_GLOB ) rt_init_context = toc(t_init_context)
        if( L_BENCH_GLOB ) t_trail_frac = tic()
        call determine_trailing_update_fraction()
        if( L_BENCH_GLOB ) rt_trail_frac = toc(t_trail_frac)
        rt_gridcorr = 0.  ! deapodization lives in the gridding backend (pair/sum restoration)
        call determine_dropped_states()
        do state = 1, params%nstates
            if( l_state_dropped(state) )then
                call carry_forward_dropped_state()
                cycle
            endif
            call assemble_state()
        enddo
        call collect_restore_timings()
        if( L_BENCH_GLOB ) t_upd_proj = tic()
        call update_project_resolution_metadata()
        if( L_BENCH_GLOB ) rt_upd_proj = toc(t_upd_proj)
        if( L_BENCH_GLOB ) t_cleanup = tic()
        call cleanup_context()
        if( L_BENCH_GLOB ) rt_cleanup = toc(t_cleanup)
        call write_benchmark()

    contains

        subroutine initialize_bench_timers()
            rt_reduce_partials           = 0.
            rt_sum_eos                   = 0.
            rt_restore_eos_and_write_fsc = 0.
            rt_restore_merged_volume     = 0.
            rt_trail_blend_accums        = 0.
            rt_trail_restored_halves     = 0.
            rt_nu_envmask                = 0.
            rt_nonuniform_filter         = 0.
            rt_init_context              = 0.
            rt_trail_frac                = 0.
            rt_gridcorr                  = 0.
            rt_upd_proj                  = 0.
            rt_cleanup                   = 0.
            rt_tot                       = 0.
            if( L_BENCH_GLOB ) t_tot = tic()
        end subroutine initialize_bench_timers

        subroutine initialize_context()
            call build%init_params_and_build_general_tbox(cline, params)
            ! the retired build_rec_eo_tbox back-filled missing per-particle
            ! projection-direction labels before assembly; preserve that contract
            if( .not. build%spproj_field%isthere('proj') )then
                call build%spproj_field%set_projs(build%eulspace)
            endif
            call set_nu_filter_report(params%part == 1)
            call even_rec%new_accumulator(params, build%spproj, expand=.false.)
            call odd_rec%new_accumulator(params, build%spproj, expand=.false.)
            call read_even_rec%new_accumulator(params, build%spproj, expand=.false.)
            call read_odd_rec%new_accumulator(params, build%spproj, expand=.false.)
            call sum_rec%new_accumulator(params, build%spproj, expand=.false.)
            numlen_part       = max(1, params%numlen)
            l_nonuniform_mode = params%l_nonuniform
            allocate(res0143s(params%nstates))
            res0143s = 0.
            allocate(nu_align_lps(params%nstates))
            nu_align_lps = 0.
            allocate(update_frac_trail_recs(params%nstates))
            update_frac_trail_recs = 1.0
            allocate(realized_update_fracs(params%nstates))
            realized_update_fracs = 1.0
            allocate(state_pops(params%nstates))
            state_pops = 0
        end subroutine initialize_context

        subroutine refresh_state_populations()
            integer :: istate
            if( .not. allocated(state_pops) ) return
            do istate = 1,params%nstates
                state_pops(istate) = build%spproj_field%get_pop(istate, 'state')
            enddo
        end subroutine refresh_state_populations

        subroutine determine_dropped_states()
            ! Multi-state reconstruction can leave a state with zero particles
            ! (e.g. independent mode where prob search reassigns all of a
            ! state's particles). Such a state produces no partial
            ! reconstructions, so re-assembling it would overwrite the previous,
            ! still-valid volume with a blank one and compute a degenerate FSC.
            ! Detect these states from the actual on-disk partial
            ! reconstructions and drop them: skip re-assembly and carry the
            ! previous volume forward as the reference.
            logical :: l_has_partials(params%nstates)
            integer :: istate
            allocate(l_state_dropped(params%nstates), source=.false.)
            if( params%nstates <= 1 ) return
            do istate = 1, params%nstates
                l_has_partials(istate) = state_has_partials(istate)
            enddo
            ! Only apply dropping in a normal iteration where at least one state
            ! was reconstructed; never drop every state.
            if( .not. any(l_has_partials) ) return
            do istate = 1, params%nstates
                if( .not. l_has_partials(istate) )then
                    l_state_dropped(istate) = .true.
                    write(logfhandle,'(A,I0,A)') &
                        &'>>> VOLASSEMBLE: STATE ', istate, &
                        &' HAS NO RECONSTRUCTED PARTICLES; DROPPING AND CARRYING PREVIOUS VOLUME FORWARD'
                endif
            enddo
        end subroutine determine_dropped_states

        logical function state_has_partials( istate ) result( l_has )
            integer, intent(in) :: istate
            type(string) :: fbody
            integer      :: part
            l_has = .false.
            do part = 1, params%nparts
                fbody = refine3D_partial_rec_fbody(istate, part, numlen_part)
                if( file_exists(fbody//'_even'//MRC_EXT) .or. file_exists(fbody//'_odd'//MRC_EXT) )then
                    l_has = .true.
                    call fbody%kill
                    return
                endif
                call fbody%kill
            enddo
        end function state_has_partials

        subroutine carry_forward_dropped_state()
            ! Preserve the existing on-disk volume, half-maps and FSC for a
            ! dropped state so downstream reference preparation still resolves a
            ! complete volume source. Nothing is written or overwritten here.
            params%vols(state)      = refine3D_state_vol_fname(state)
            params%vols_even(state) = refine3D_state_halfvol_fname(state, 'even')
            params%vols_odd(state)  = refine3D_state_halfvol_fname(state, 'odd')
            res0143s(state)         = 0.
        end subroutine carry_forward_dropped_state

        !> The realized fractions f (what the current partials actually contain)
        !! are always computed; the applied fractions u (the requested map-update
        !! weight) default to f and are overridden by ufrac_trec for single-state
        !! runs. The accumulator blend needs both: u sets the previous-chain
        !! decay and f normalizes the current partials to weight u.
        subroutine determine_trailing_update_fraction()
            update_frac_trail_recs = 1.0
            realized_update_fracs  = 1.0
            if( .not. params%l_trail_rec ) return
            call build%spproj%read_segment(params%oritype, params%projfile)
            call build%spproj%os_ptcl3D%get_state_update_fracs(params%nstates, realized_update_fracs)
            update_frac_trail_recs = realized_update_fracs
            if( params%l_ufrac_trec_defined )then
                if( params%nstates == 1 )then
                    update_frac_trail_recs(1) = params%ufrac_trec
                    write(logfhandle,'(A,F8.4)') &
                        '>>> VOLASSEMBLE: USING EXPLICIT SINGLE-STATE UFRAC_TREC = ', params%ufrac_trec
                else
                    THROW_WARN('ufrac_trec is ignored for multi-state trailing reconstruction; using realized per-state update fractions')
                endif
            endif
        end subroutine determine_trailing_update_fraction

        subroutine assemble_state()
            call restore_state_from_parts(params, build, cline, even_rec, odd_rec, read_even_rec, read_odd_rec, &
                &sum_rec, state, numlen_part, &
                &update_frac_trail_recs(state), realized_update_fracs(state), &
                &vol_prev_even, vol_prev_odd, vol_merged, &
                &vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd, &
                &volname, eonames, res05, res0143s(state), restore_timings)
            params%vols(state)      = volname
            params%vols_even(state) = eonames(1)
            params%vols_odd(state)  = eonames(2)
            call postprocess_state()
            call volname%kill
        end subroutine assemble_state

        subroutine postprocess_state()
            integer :: which_iter
            which_iter = 1
            if( cline%defined('which_iter') ) which_iter = params%which_iter
            call plan_state_postprocess(params, state, which_iter, pp_plan)
            if( pp_plan%l_nu_envmask_incompatible )then
                write(logfhandle,'(A,1X,A)') &
                    &'>>> Existing NU evidence envelope incompatible with current box/sampling, regenerating:', &
                    &pp_plan%nu_envmask_file%to_char()
            endif
            if( l_nonuniform_mode ) call run_state_nonuniform_filter()
        end subroutine postprocess_state

        subroutine run_state_nonuniform_filter()
            if( L_BENCH_GLOB ) t_nonuniform_filter = tic()
            call setup_nonuniform_filter()
            call release_nonuniform_aux_inputs()
            call optimize_nu_cutoff_finds()
            call refine_nonuniform_filter_bank()
            call generate_state_nu_envmask()
            call release_nonuniform_base_inputs()
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
            call log_nonuniform_filter_stats()
            call write_nonuniform_outputs()
            call record_nu_alignment_lowpass_limit()
            call cleanup_nonuniform_state()
            if( L_BENCH_GLOB ) rt_nonuniform_filter = rt_nonuniform_filter + toc(t_nonuniform_filter)
        end subroutine run_state_nonuniform_filter

        subroutine generate_state_nu_envmask()
            type(nu_envmask_params) :: envp
            type(nu_envmask_stats)  :: envstats
            logical, allocatable    :: l_env(:,:,:)
            integer                 :: n_ccs, n_ccs_kept, grow_px, edge_px
            if( pp_plan%nu_envmask_action /= NU_ENVMASK_ACTION_REGENERATE ) return
            if( L_BENCH_GLOB ) t_nu_envmask = tic()
            envp%nsigma      = params%nu_msk_sig
            envp%beta        = NU_ENVMASK_BETA
            envp%dens_weight = NU_ENVMASK_DENS_WEIGHT
            envp%lp_smooth   = params%amsklp
            envp%l_relative  = NU_ENVMASK_RELATIVE
            call nu_evidence_envelope(envp, l_env, envstats)
            call print_nu_envmask_stats(envstats)
            grow_px = max(1, nint(NU_ENVMASK_GROW_A / vol_nu_base_even%get_smpd()))
            edge_px = max(1, nint(NU_ENVMASK_EDGE_A / vol_nu_base_even%get_smpd()))
            call nu_envmask%kill_bimg
            call nu_envmask%envmask3D_from_lmask(l_env, vol_nu_base_even%get_smpd(), grow_px, edge_px, &
                &NU_ENVMASK_MINVOL_FRAC, .true., n_ccs, n_ccs_kept)
            ! an empty evidence field returns without constructing the mask
            ! image; writing it would abort on invalid MRC dimensions. Skip the
            ! write -- every envelope consumer handles absence.
            if( n_ccs_kept < 1 )then
                THROW_WARN('NU evidence envelope is empty; no envelope mask written this iteration')
            else
                call nu_envmask%write(pp_plan%nu_envmask_file, del_if_exists=.true.)
                call wait_for_closure(pp_plan%nu_envmask_file)
                write(logfhandle,'(A,I0,A,1X,A)') &
                    &'>>> NU EVIDENCE ENVELOPE: STATE ', state, ', MASK', pp_plan%nu_envmask_file%to_char()
            endif
            ! One greppable line per state per cycle. The envelope is allowed to
            ! shrink as resolution improves, but reference masking suppresses
            ! whatever it excludes, so a monotonically falling occupancy is the
            ! signal that it is clipping rather than tightening.
            write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') &
                &'>>> NU ENVELOPE OCCUPANCY: STATE ', state, ', SUPPORT FRACTION ', &
                &envstats%pct_signal, ' %, COMPONENTS KEPT ', n_ccs_kept, ' OF ', n_ccs
            if( allocated(l_env) ) deallocate(l_env)
            if( L_BENCH_GLOB ) rt_nu_envmask = rt_nu_envmask + toc(t_nu_envmask)
        end subroutine generate_state_nu_envmask

        subroutine setup_nonuniform_filter()
            integer :: n_highres_steps
            real    :: aux_resolution
            n_highres_steps = nu_highres_steps_for_state()
            call cleanup_nu_aux_images()
            if( use_static_nu_aux_replacement() )then
                allocate(nu_aux_even(1), nu_aux_odd(1))
                call nu_aux_even(1)%copy(vol_nu_aux_even)
                call nu_aux_odd(1)%copy(vol_nu_aux_odd)
                aux_resolution = nu_aux_effective_resolution()
                call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, params%mskdiam, [aux_resolution], &
                    &nu_aux_even, nu_aux_odd, n_highres_steps=n_highres_steps)
            else
                call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, params%mskdiam, [real ::], &
                    &n_highres_steps=n_highres_steps)
            endif
        end subroutine setup_nonuniform_filter

        logical function use_static_nu_aux_replacement() result(l_use_aux)
            l_use_aux = params%l_ml_reg .and. .not. params%l_nu_refine
        end function use_static_nu_aux_replacement

        real function nu_aux_effective_resolution() result(aux_resolution)
            aux_resolution = res0143s(state)
            if( params%l_lpset .and. params%lp > TINY )then
                if( NU_DEV_OUTPUT .and. params%part == 1 .and. aux_resolution > params%lp + TINY )then
                    write(logfhandle,'(A,F8.3,A,F8.3,A)') &
                        &'>>> NU auxiliary effective resolution clamped by matching low-pass: FSC ', &
                        &aux_resolution, ' A; matching LP ', params%lp, ' A'
                endif
                aux_resolution = min(aux_resolution, params%lp)
            endif
        end function nu_aux_effective_resolution

        subroutine refine_nonuniform_filter_bank()
            type(nu_highres_extension_stats) :: ext_stats
            integer :: n_highres_steps, n_accepted_this_iteration
            if( .not. params%l_nu_refine ) return
            n_highres_steps = nu_highres_steps_for_state()
            n_accepted_this_iteration = 0
            do
                call extend_nu_filter_highres_shell_next(vol_nu_base_even, vol_nu_base_odd, stats=ext_stats)
                if( .not. ext_stats%attempted )then
                    if( NU_DEV_OUTPUT .and. params%part == 1 )then
                        if( ext_stats%n_mask == 0 )then
                            write(logfhandle,'(A)') &
                                &'>>> NU high-resolution extension stopped: empty NU refinement mask'
                        else if( ext_stats%n_tested == 0 )then
                            write(logfhandle,'(A,F8.3,A,I0,A)') &
                                &'>>> NU high-resolution extension stopped: no frontier voxels at current finest label ', &
                                &ext_stats%old_limit, ' A (k=', ext_stats%old_find, ')'
                        else
                            write(logfhandle,'(A,F8.3,A,I0,A)') &
                                &'>>> NU high-resolution extension stopped: no valid next shell after ', &
                                &ext_stats%old_limit, ' A (k=', ext_stats%old_find, ')'
                        endif
                    endif
                    exit
                endif
                if( .not. ext_stats%applied      ) exit
                if( .not. ext_stats%promote_next ) exit
                n_accepted_this_iteration = n_accepted_this_iteration + 1
            end do
            if( n_accepted_this_iteration > 0 )then
                call refine_nu_extension_filtmap_ordered_labels
                n_highres_steps = get_nu_filtmap_highres_shell_depth()
                call write_nu_highres_steps_for_state(n_highres_steps)
                if( NU_DEV_OUTPUT .and. params%part == 1 )then
                    write(logfhandle,'(A,I0,A,I0)') &
                        &'>>> NU high-resolution extension accepted shell steps this iteration: ', &
                        &n_accepted_this_iteration, '; promoted depth for next iteration: ', n_highres_steps
                endif
            endif
        end subroutine refine_nonuniform_filter_bank

        integer function nu_highres_steps_for_state() result(nsteps)
            type(string) :: fname
            integer :: funit, io_stat
            nsteps = 0
            if( .not. params%l_nu_refine ) return
            if( params%startit <= 1 .and. params%which_iter <= params%startit )then
                call write_nu_highres_steps_for_state(0)
                return
            endif
            fname = nu_highres_steps_fname()
            if( .not.file_exists(fname) )then
                call fname%kill
                return
            endif
            call fopen(funit, status='OLD', action='READ', file=fname, iostat=io_stat)
            if( io_stat == 0 )then
                read(funit, *, iostat=io_stat) nsteps
                call fclose(funit)
            endif
            if( io_stat /= 0 ) nsteps = 0
            nsteps = max(0, nsteps)
            call fname%kill
        end function nu_highres_steps_for_state

        subroutine write_nu_highres_steps_for_state(nsteps)
            integer, intent(in) :: nsteps
            type(string) :: fname
            integer :: funit, io_stat
            if( .not. params%l_nu_refine ) return
            fname = nu_highres_steps_fname()
            call fopen(funit, status='REPLACE', action='WRITE', file=fname, iostat=io_stat)
            if( io_stat == 0 )then
                write(funit,'(I0)') max(0, nsteps)
                call fclose(funit)
            else
                write(logfhandle,'(A,1X,A)') '>>> WARNING: failed to write NU high-resolution depth file:', &
                    &fname%to_char()
            endif
            call fname%kill
        end subroutine write_nu_highres_steps_for_state

        function nu_highres_steps_fname() result(fname)
            type(string) :: fname
            fname = 'nu_highres_depth_state'//int2str_pad(state,2)//'.txt'
        end function nu_highres_steps_fname

        subroutine log_nonuniform_filter_stats()
            call print_nu_filtmap_lowpass_stats()
            if( NU_DEV_OUTPUT .and. params%part == 1 ) call analyze_filtmap_neighbor_continuity()
        end subroutine log_nonuniform_filter_stats

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
            real :: selected_lp
            if( .not. params%l_nonuniform ) return
            selected_lp = get_nu_filtmap_finest_selected_lp()
            if( selected_lp <= TINY ) return
            nu_align_lps(state) = selected_lp
            if( NU_DEV_OUTPUT .and. params%part == 1 )then
                write(logfhandle,'(A,I0,A,F8.3,A)') &
                    &'>>> NU filter state ', state, ' matching low-pass limit for next iteration: ', selected_lp, ' A'
            endif
        end subroutine record_nu_alignment_lowpass_limit

        subroutine cleanup_nonuniform_state()
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
            call vol_nu_aux_even%kill
            call vol_nu_aux_odd%kill
            call cleanup_nu_aux_images()
            call cleanup_nu_filter()
        end subroutine cleanup_nonuniform_state

        subroutine release_nonuniform_aux_inputs()
            call cleanup_nu_aux_images()
            call vol_nu_aux_even%kill
            call vol_nu_aux_odd%kill
        end subroutine release_nonuniform_aux_inputs

        subroutine release_nonuniform_base_inputs()
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
        end subroutine release_nonuniform_base_inputs

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

        subroutine collect_restore_timings()
            rt_reduce_partials           = restore_timings%reduce_partials
            rt_sum_eos                   = restore_timings%sum_eos
            rt_restore_eos_and_write_fsc = restore_timings%restore_eos_and_write_fsc
            rt_restore_merged_volume     = restore_timings%restore_merged_volume
            rt_trail_blend_accums        = restore_timings%trail_blend_accums
            rt_trail_restored_halves     = restore_timings%trail_restored_halves
        end subroutine collect_restore_timings

        subroutine update_project_resolution_metadata()
            integer :: iptcl, istate
            call refresh_state_populations()
            if( params%nstates == 1 )then
                call build%spproj_field%set_all2single('res', res0143s(1))
            else
                do iptcl = 1, build%spproj_field%get_noris()
                    istate = build%spproj_field%get_state(iptcl)
                    if( istate > 0 .and. istate <= params%nstates )then
                        call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                    endif
                enddo
            endif
            call update_project_nu_alignment_lowpass()
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        end subroutine update_project_resolution_metadata

        subroutine update_project_nu_alignment_lowpass()
            logical :: l_included(params%nstates)
            real    :: align_lp
            integer :: istate, selected_state
            if( .not. l_nonuniform_mode ) return
            if( .not. allocated(nu_align_lps) ) return
            if( .not. allocated(state_pops) ) return
            ! Match the classical multi-state policy: the best resolved
            ! populated state determines the single global matching bandwidth.
            l_included = (nu_align_lps > TINY) .and. (state_pops > 0)
            if( NU_DEV_OUTPUT .and. params%part == 1 .and. params%nstates > 1 ) call log_nu_alignment_lowpass_summary(l_included)
            if( .not. any(l_included) )then
                if( any(nu_align_lps > TINY) )then
                    write(logfhandle,'(A)') &
                        &'>>> WARNING: no populated state has a valid NU filter matching low-pass limit'
                endif
                return
            endif
            align_lp       = minval(nu_align_lps, mask=l_included)
            selected_state = 0
            do istate = 1,params%nstates
                if( l_included(istate) .and. abs(nu_align_lps(istate) - align_lp) <= TINY )then
                    selected_state = istate
                    exit
                endif
            enddo
            call build%spproj_field%set_all2single('lp', align_lp)
            if( NU_DEV_OUTPUT .and. params%part == 1 )then
                write(logfhandle,'(A,I0,A,F8.3,A)') &
                    &'>>> NU filter project matching low-pass limit from state ', selected_state, ': ', align_lp, ' A'
            endif
        end subroutine update_project_nu_alignment_lowpass

        subroutine log_nu_alignment_lowpass_summary(l_included)
            logical, intent(in) :: l_included(:)
            integer :: istate
            write(logfhandle,'(A)') '>>> NU filter multi-state matching low-pass candidates'
            write(logfhandle,'(A)') '    State       Pop   FSC(A)   NU LP(A)   Used'
            do istate = 1,params%nstates
                write(logfhandle,'(I9,I10,F9.3,F10.3,5X,A)') &
                    &istate, state_pops(istate), res0143s(istate), nu_align_lps(istate), &
                    &merge('yes', 'no ', l_included(istate))
            enddo
        end subroutine log_nu_alignment_lowpass_summary

        subroutine cleanup_context()
            call build%kill_general_tbox
            call even_rec%kill
            call odd_rec%kill
            call read_even_rec%kill
            call read_odd_rec%kill
            call sum_rec%kill
            call vol_prev_even%kill
            call vol_prev_odd%kill
            call vol_merged%kill
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
            call vol_nu_aux_even%kill
            call vol_nu_aux_odd%kill
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call cleanup_nu_aux_images()
            if( allocated(state_pops)             ) deallocate(state_pops)
            if( allocated(l_state_dropped)        ) deallocate(l_state_dropped)
            if( allocated(res0143s)               ) deallocate(res0143s)
            if( allocated(nu_align_lps)           ) deallocate(nu_align_lps)
            if( allocated(update_frac_trail_recs) ) deallocate(update_frac_trail_recs)
            if( allocated(realized_update_fracs)  ) deallocate(realized_update_fracs)
            call cleanup_nu_filter()
            call nu_envmask%kill_bimg
            call pp_plan%nu_envmask_file%kill
            call volname%kill
            call eonames(1)%kill
            call eonames(2)%kill
            call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
            call simple_touch('VOLASSEMBLE_FINISHED')
        end subroutine cleanup_context

        subroutine write_benchmark()
            type(string) :: benchfname
            integer :: fnr
            if( .not. L_BENCH_GLOB ) return
            rt_tot         = toc(t_tot)
            benchfname = refine3D_volassemble_bench_fname(params%which_iter)
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a)')    'volassemble assembly mode             : volume'
            write(fnr,'(a,i0)') 'volassemble nspace                    : ', params%nspace
            write(fnr,'(a,i0)') 'volassemble nstates                   : ', params%nstates
            write(fnr,'(a,i0)') 'volassemble kfrom                     : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'volassemble kto                       : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f0.2)') 'volassemble reduce_partials           :', rt_reduce_partials
            write(fnr,'(a,1x,f0.2)') 'volassemble sum_eos                   :', rt_sum_eos
            write(fnr,'(a,1x,f0.2)') 'volassemble restore_eos_and_write_fsc :', &
                rt_restore_eos_and_write_fsc
            write(fnr,'(a,1x,f0.2)') 'volassemble restore_merged_volume     :', rt_restore_merged_volume
            write(fnr,'(a,1x,f0.2)') 'volassemble trail_blend_accums        :', rt_trail_blend_accums
            write(fnr,'(a,1x,f0.2)') 'volassemble trail_restored_halves     :', rt_trail_restored_halves
            write(fnr,'(a,1x,f0.2)') 'volassemble nu_evidence_envelope      :', rt_nu_envmask
            write(fnr,'(a,1x,f0.2)') 'volassemble nonuniform_filter         :', rt_nonuniform_filter
            write(fnr,'(a,1x,f0.2)') 'volassemble init_context              :', rt_init_context
            write(fnr,'(a,1x,f0.2)') 'volassemble trail_frac_read           :', rt_trail_frac
            write(fnr,'(a,1x,f0.2)') 'volassemble gridcorr_prep             :', rt_gridcorr
            write(fnr,'(a,1x,f0.2)') 'volassemble update_proj_metadata      :', rt_upd_proj
            write(fnr,'(a,1x,f0.2)') 'volassemble cleanup                   :', rt_cleanup
            write(fnr,'(a,1x,f0.2)') 'volassemble total time                :', rt_tot
            write(fnr,'(a,1x,f0.2)') 'volassemble % accounted for           :', &
                &((rt_reduce_partials + rt_sum_eos + rt_restore_eos_and_write_fsc +               &
                &  rt_restore_merged_volume + rt_trail_blend_accums + rt_trail_restored_halves +  &
                &  rt_nonuniform_filter +                                                         &
                &  rt_init_context + rt_trail_frac + rt_gridcorr + rt_upd_proj + rt_cleanup)      &
                & / rt_tot) * 100.
            call fclose(fnr)
            call benchfname%kill
        end subroutine write_benchmark

    end subroutine exec_volassemble

end module simple_commanders_rec_distr
