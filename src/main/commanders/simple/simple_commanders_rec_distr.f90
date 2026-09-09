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
        use simple_reconstructor, only: reconstructor, gridding_half_restore
        use simple_halfmap_diagnostics, only: halfmap_diagnostics_result, evaluate_halfmap_pair, &
            &write_halfmap_diagnostics, write_support_provenance
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
        type(halfmap_diagnostics_result) :: pair_diagnostics
        real    :: weight_prev
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
        ! The shipped halves and merged volume carry the soft spherical support
        ! at msk_crop (restore_gridding_pair, restore_merged_volume; the legacy
        ! trailing blend mixes two such volumes). Record it beside the volume
        ! so postprocess does not mask again; the PCG base warm-start selector
        ! ignores the gridding kind.
        call write_support_provenance(volname, .false., 'gridding')
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
            type(fsc_area_score_result)      :: cones_fsc
            type(halfmap_diagnostics_result) :: prev_diagnostics
            if( L_BENCH_GLOB ) t_restore_eos = tic()
            if( params%l_trail_rec .and. .not. l_trail_chain )then
                ! Bootstrap iteration: the chain has no information yet, so the
                ! previous halfmaps provide the FSC prior for regularization,
                ! exactly as the legacy volume-domain trailing did. The
                ! previous half maps are final deapodized real-space volumes;
                ! they satisfy the common evaluator's real-space contract
                ! directly, with no representation adapter.
                call read_previous_halfmaps()
                if( params%conical_fsc == 'yes' )then
                    call calc_gridding_pair_diagnostics(params, vol_prev_even, vol_prev_odd, &
                        &state, prev_diagnostics, cones=cones_fsc)
                    call restore_gridding_pair(params, even_rec, odd_rec, state, &
                        &eonames(1), eonames(2), pair_diagnostics, fsc_in=prev_diagnostics%fsc, &
                        &cfar_in=prev_diagnostics%cfar, cones_in=cones_fsc)
                    call cones_fsc%kill
                else
                    call calc_gridding_pair_diagnostics(params, vol_prev_even, vol_prev_odd, &
                        &state, prev_diagnostics)
                    call restore_gridding_pair(params, even_rec, odd_rec, state, &
                        &eonames(1), eonames(2), pair_diagnostics, fsc_in=prev_diagnostics%fsc, &
                        &cfar_in=prev_diagnostics%cfar)
                endif
                call prev_diagnostics%kill
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
            call write_halfmap_diagnostics(pair_diagnostics, params%box_crop, params%smpd_crop, fsc_txt_file)
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
            ! the same soft spherical support as the halves and the PCG solve
            call build%vol%mask3D_soft(params%msk_crop, backgr=0.)
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
        end subroutine cleanup_restore_state

        ! GRIDDING PAIR ACCUMULATOR I/O AND RESTORATION
        ! These procedures coordinate workflow artifacts and type-bound domain
        ! operations for this host only; the numerical kernels stay type-bound
        ! on reconstructor and the diagnostics in simple_halfmap_diagnostics.
        ! Arguments stay explicit so mutation and lifecycle remain visible.

        !> Read an explicit even/odd gridding-accumulator artifact. Previous
        !! smaller grids retain the legacy update-fraction zero-padding behavior.
        !! A missing artifact resets the pair (legacy partial semantics) unless the
        !! caller marks the artifact required, in which case it is a hard error.
        subroutine read_gridding_pair_accumulators( params, even_rec, odd_rec, fbody, required )
            use simple_imgfile, only: imgfile
            class(parameters),    intent(in)    :: params
            class(reconstructor), intent(inout) :: even_rec, odd_rec
            class(string),        intent(in)    :: fbody
            logical, optional,    intent(in)    :: required
            type(string)      :: even_vol, even_rho, odd_vol, odd_rho
            type(image)       :: prev_vol_e, prev_vol_o
            type(imgfile)     :: ioimg_e, ioimg_o
            real, allocatable :: rho_e(:,:,:), rho_o(:,:,:)
            integer :: current_ldim(3), cshape(3), prev_ldim(3)
            integer :: fhandle_rho_e, fhandle_rho_o, i, ierr, dummy
            real    :: current_smpd, prev_smpd
            logical :: here(4), l_pad_with_zeros
            current_ldim = even_rec%get_ldim()
            current_smpd = even_rec%get_smpd()
            if( any(odd_rec%get_ldim() /= current_ldim) )then
                THROW_HARD('gridding pair accumulator dimensions do not match')
            endif
            if( abs(odd_rec%get_smpd() - current_smpd) > TINY )then
                THROW_HARD('gridding pair accumulator sampling does not match')
            endif
            even_vol = fbody//'_even'//MRC_EXT
            even_rho = string('rho_')//fbody//'_even'//MRC_EXT
            odd_vol  = fbody//'_odd'//MRC_EXT
            odd_rho  = string('rho_')//fbody//'_odd'//MRC_EXT
            here(1)  = file_exists(even_vol)
            here(2)  = file_exists(even_rho)
            here(3)  = file_exists(odd_vol)
            here(4)  = file_exists(odd_rho)
            if( all(here) )then
                l_pad_with_zeros = .false.
                if( params%l_update_frac )then
                    call find_ldim_nptcls(even_vol, prev_ldim, dummy)
                    prev_smpd = current_smpd
                    if( prev_ldim(1) == current_ldim(1) )then
                        ! matching grids require no compatibility transform
                    elseif( prev_ldim(1) > current_ldim(1) )then
                        THROW_HARD('previous gridding accumulator is larger than the current grid')
                    else
                        l_pad_with_zeros = .true.
                    endif
                endif
                call fopen(fhandle_rho_e, file=even_rho, status='OLD', action='READ', access='STREAM', iostat=ierr)
                call fileiochk('read_gridding_pair_accumulators opening '//even_rho%to_char(), ierr)
                call fopen(fhandle_rho_o, file=odd_rho, status='OLD', action='READ', access='STREAM', iostat=ierr)
                call fileiochk('read_gridding_pair_accumulators opening '//odd_rho%to_char(), ierr)
                if( l_pad_with_zeros )then
                    call ioimg_e%open(even_vol, prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                    call ioimg_o%open(odd_vol, prev_ldim, prev_smpd, formatchar='M', readhead=.false., rwaction='read')
                    call prev_vol_e%new(prev_ldim, prev_smpd)
                    call prev_vol_o%new(prev_ldim, prev_smpd)
                    cshape = [fdim(prev_ldim(1)), prev_ldim(2), prev_ldim(3)]
                    allocate(rho_e(1:cshape(1),1:cshape(2),1:cshape(3)), &
                        &rho_o(1:cshape(1),1:cshape(2),1:cshape(3)))
                    call even_rec%reset
                    call odd_rec%reset
                    !$omp parallel do default(shared) private(i,ierr) schedule(static) num_threads(4)
                    do i = 1, 4
                        select case(i)
                            case(1)
                                call prev_vol_e%read_raw_mrc(ioimg_e)
                            case(2)
                                call prev_vol_o%read_raw_mrc(ioimg_o)
                            case(3)
                                read(fhandle_rho_e, pos=1, iostat=ierr) rho_e
                                if( ierr /= 0 ) call fileiochk('read_gridding_pair_accumulators reading even rho', ierr)
                            case(4)
                                read(fhandle_rho_o, pos=1, iostat=ierr) rho_o
                                if( ierr /= 0 ) call fileiochk('read_gridding_pair_accumulators reading odd rho', ierr)
                        end select
                    enddo
                    !$omp end parallel do
                    call even_rec%pad_with_zeros(prev_vol_e, rho_e)
                    call odd_rec%pad_with_zeros(prev_vol_o, rho_o)
                    call prev_vol_e%kill
                    call prev_vol_o%kill
                    deallocate(rho_e, rho_o)
                else
                    call ioimg_e%open(even_vol, current_ldim, current_smpd, formatchar='M', &
                        &readhead=.false., rwaction='read')
                    call ioimg_o%open(odd_vol, current_ldim, current_smpd, formatchar='M', &
                        &readhead=.false., rwaction='read')
                    !$omp parallel do default(shared) private(i) schedule(static) num_threads(4)
                    do i = 1, 4
                        select case(i)
                            case(1)
                                call even_rec%read_raw_mrc(ioimg_e)
                            case(2)
                                call odd_rec%read_raw_mrc(ioimg_o)
                            case(3)
                                call even_rec%read_raw_rho(fhandle_rho_e)
                            case(4)
                                call odd_rec%read_raw_rho(fhandle_rho_o)
                        end select
                    enddo
                    !$omp end parallel do
                endif
                call ioimg_e%close
                call ioimg_o%close
                call fclose(fhandle_rho_e)
                call fclose(fhandle_rho_o)
            else
                if( present(required) )then
                    if( required )then
                        THROW_HARD('required gridding pair accumulator artifact is incomplete: '//fbody%to_char())
                    endif
                endif
                call even_rec%reset
                call odd_rec%reset
            endif
        end subroutine read_gridding_pair_accumulators

        !> Write an explicit even/odd gridding-accumulator artifact using the
        !! established numerator/rho filenames.
        subroutine write_gridding_pair_accumulators( even_rec, odd_rec, fbody )
            class(reconstructor), intent(inout) :: even_rec, odd_rec
            class(string),        intent(in)    :: fbody
            call even_rec%write_raw_accum(fbody//'_even'//MRC_EXT, &
                &string('rho_')//fbody//'_even'//MRC_EXT)
            call odd_rec%write_raw_accum(fbody//'_odd'//MRC_EXT, &
                &string('rho_')//fbody//'_odd'//MRC_EXT)
        end subroutine write_gridding_pair_accumulators

        !> Restore one explicit gridding half pair. This owns backend-specific
        !! base/replay mechanics but no composite lifetime, partial reduction, or
        !! trailing-chain policy.
        subroutine restore_gridding_pair( params, even_rec, odd_rec, state, fname_even, fname_odd, &
            &diagnostics, fsc_in, cfar_in, cones_in )
            use simple_fsc, only: fsc_area_score_result
            class(parameters),                      intent(in)    :: params
            class(reconstructor),                   intent(inout) :: even_rec, odd_rec
            integer,                                intent(in)    :: state
            class(string),                          intent(in)    :: fname_even, fname_odd
            type(halfmap_diagnostics_result),       intent(out)   :: diagnostics
            real, optional,                         intent(in)    :: fsc_in(:)
            real, optional,                         intent(in)    :: cfar_in
            class(fsc_area_score_result), optional, intent(inout) :: cones_in
            type(gridding_half_restore) :: even_restore, odd_restore
            type(fsc_area_score_result) :: cones_fsc
            real,     allocatable :: res(:)
            real                  :: smpd, fny
            integer               :: box, filtsz
            logical               :: l_have_fsc
            box    = params%box_crop
            smpd   = params%smpd_crop
            filtsz = fdim(box) - 1
            fny    = 2. * smpd
            res    = get_resarr(box, smpd)
            if( present(fsc_in) )then
                if( .not. present(cfar_in) ) THROW_HARD('cfar_in must accompany an input gridding FSC')
                if( size(fsc_in) /= filtsz ) THROW_HARD('input FSC size does not match gridding reconstruction')
                allocate(diagnostics%fsc(filtsz), source=fsc_in)
                diagnostics%cfar = cfar_in
                l_have_fsc = .true.
                if( params%l_ml_reg .and. (params%conical_fsc == 'yes') )then
                    if( .not. present(cones_in) )then
                        THROW_HARD('cones_in must be provided if conical regularization is enabled')
                    endif
                endif
            else
                allocate(diagnostics%fsc(filtsz), source=0.)
                l_have_fsc = .false.
            endif
            call even_restore%new(even_rec)
            call odd_restore%new(odd_rec)
            ! Every gridding product is deapodized and then given the soft
            ! spherical support the PCG solve installs (mask3D_soft at
            ! msk_crop), so the two backends' halves are equivalent and the FSC
            ! is computed on the halves as shipped, with no mask of its own
            ! (2026-09-09). The inverse envelope's rim gain, which the legacy
            ! undeapodized FSC avoided, is zeroed by the support.
            if( params%l_ml_reg )then
                call even_rec%restore_base(even_restore%base, preserve_numerator=.true.)
                call even_restore%finalize_from_base(even_rec)
                call even_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call even_restore%final%write(add2fbody(fname_even,MRC_EXT,'_unfil'), del_if_exists=.true.)
                call odd_rec%restore_base(odd_restore%base, preserve_numerator=.true.)
                call odd_restore%finalize_from_base(odd_rec)
                call odd_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call odd_restore%final%write(add2fbody(fname_odd,MRC_EXT,'_unfil'), del_if_exists=.true.)
                if( .not. l_have_fsc )then
                    call calc_gridding_pair_diagnostics(params, even_restore%final, odd_restore%final, &
                        &state, diagnostics, cones=cones_fsc)
                endif
                call even_restore%final%kill
                call odd_restore%final%kill
                ! Regularization
                if( l_have_fsc )then
                    if( params%conical_fsc == 'yes' )then
                        call even_rec%add_conical_invtausq2rho(cones_in)
                        call odd_rec%add_conical_invtausq2rho(cones_in)
                    else
                        call even_rec%add_invtausq2rho(diagnostics%fsc)
                        call odd_rec%add_invtausq2rho(diagnostics%fsc)
                    endif
                else
                    if( params%conical_fsc == 'yes' )then
                        call even_rec%add_conical_invtausq2rho(cones_fsc)
                        call odd_rec%add_conical_invtausq2rho(cones_fsc)
                    else
                        call even_rec%add_invtausq2rho(diagnostics%fsc)
                        call odd_rec%add_invtausq2rho(diagnostics%fsc)
                    endif
                endif
                ! regularized halves: density correction, deapodization, support, write
                call even_restore%base%kill
                call even_restore%prepare_final(even_rec)
                call even_rec%restore_final(even_restore%final, preserve_numerator=.true.)
                call even_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call even_restore%final%write(fname_even, del_if_exists=.true.)
                call even_restore%final%kill
                call odd_restore%base%kill
                call odd_restore%prepare_final(odd_rec)
                call odd_rec%restore_final(odd_restore%final, preserve_numerator=.true.)
                call odd_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call odd_restore%final%write(fname_odd, del_if_exists=.true.)
                call odd_restore%final%kill
            else
                call even_rec%restore_base(even_restore%base)
                call odd_rec%restore_base(odd_restore%base)
                call even_restore%finalize_from_base(even_rec)
                call even_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call even_restore%final%write(fname_even, del_if_exists=.true.)
                call odd_restore%finalize_from_base(odd_rec)
                call odd_restore%final%mask3D_soft(params%msk_crop, backgr=0.)
                call odd_restore%final%write(fname_odd, del_if_exists=.true.)
                if( .not. l_have_fsc )then
                    call calc_gridding_pair_diagnostics(params, even_restore%final, odd_restore%final, &
                        &state, diagnostics, cones=cones_fsc)
                endif
                call even_restore%final%kill
                call odd_restore%final%kill
            endif
            ! save, get & print resolution
            call arr2file(diagnostics%fsc, refine3D_fsc_fname(state))
            call get_resolution(diagnostics%fsc, res, diagnostics%res_fsc05, diagnostics%res_fsc0143)
            diagnostics%res_fsc05   = max(diagnostics%res_fsc05, fny)
            diagnostics%res_fsc0143 = max(diagnostics%res_fsc0143, fny)
            deallocate(res)
            call even_restore%kill
            call odd_restore%kill
            call cones_fsc%kill
        end subroutine restore_gridding_pair

        !> Gridding adapter for the backend-neutral half-map evaluator: builds
        !! the merged average and writes the automask artifact on the envfsc
        !! path. The evaluator applies no mask (2026-09-09): the ordinary
        !! restoration path passes the deapodized halves that already carry
        !! the soft spherical support at msk_crop (identical to the PCG solve
        !! support), and the trailing bootstrap passes the previous final half
        !! maps, which were shipped under the same contract.
        subroutine calc_gridding_pair_diagnostics( params, even, odd, state, diagnostics, cones )
            use simple_fsc, only: fsc_area_score_result
            class(parameters),                      intent(in)    :: params
            class(image),                           intent(in)    :: even, odd
            integer,                                intent(in)    :: state
            type(halfmap_diagnostics_result),       intent(out)   :: diagnostics
            class(fsc_area_score_result), optional, intent(inout) :: cones
            type(image) :: average, envmask
            call average%copy(even)
            call average%add(odd)
            call average%mul(0.5)
            if( params%l_envfsc )then
                call evaluate_halfmap_pair(params, state, even, odd, average, diagnostics, &
                    &envmask=envmask, cones=cones)
                call envmask%write(string(AUTOMASK_FBODY//int2str_pad(state,2)//MRC_EXT))
                call envmask%kill
            else
                call evaluate_halfmap_pair(params, state, even, odd, average, diagnostics, cones=cones)
            endif
            call average%kill
        end subroutine calc_gridding_pair_diagnostics

    end subroutine restore_state_from_parts

    !> PCG-path project handoff of the NU matching low-pass, the dual of
    !! update_project_nu_alignment_lowpass in the gridding volassemble: the
    !! _nu_filt products themselves are written by nonuniform_filter_state
    !! inside the PCG master (both backends run the same competition).
    subroutine filter_pcg_nonuniform_maps( params, build, l_trail_bootstrap, nu_align_lps )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        logical,          intent(in)    :: l_trail_bootstrap(:)
        real, optional,   intent(in)    :: nu_align_lps(:)
        integer, allocatable :: state_pops(:)
        logical, allocatable :: l_included(:)
        integer :: state
        real    :: align_lp
        if( .not. params%l_nonuniform ) return
        if( size(l_trail_bootstrap) /= params%nstates ) &
            &THROW_HARD('PCG nonuniform trailing-bootstrap state input has invalid size')
        if( .not. present(nu_align_lps) ) &
            &THROW_HARD('NU on the PCG backend requires the matching low-pass handoff from the NU filter')
        if( size(nu_align_lps) /= params%nstates ) &
            &THROW_HARD('NU matching low-pass handoff has invalid size')
        allocate(state_pops(params%nstates), l_included(params%nstates))
        do state = 1, params%nstates
            state_pops(state) = build%spproj_field%get_pop(state, 'state')
        enddo
        l_included = (state_pops > 0) .and. (nu_align_lps > TINY)
        if( any(l_included) )then
            align_lp = minval(nu_align_lps, mask=l_included)
            write(logfhandle,'(A,F8.3,A)') '>>> NU filter project matching low-pass limit: ', align_lp, ' A'
            call build%spproj_field%set_all2single('lp',     align_lp)
            call build%spproj_field%set_all2single('lp_est', align_lp)
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        else
            write(logfhandle,'(A)') '>>> WARNING: no populated state has a valid NU filter matching low-pass limit'
        endif
        deallocate(state_pops, l_included)
    end subroutine filter_pcg_nonuniform_maps

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor,    only: reconstructor
        use simple_nu_filter,        only: set_nu_filter_report, NU_DEV_OUTPUT
        use simple_nu_state_filter,  only: nonuniform_filter_state, nu_state_filter_timings, &
            &nu_static_aux_replacement
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters), target      :: params
        type(builder)                 :: build
        type(reconstructor)           :: even_rec, odd_rec, read_even_rec, read_odd_rec, sum_rec
        type(image)                   :: vol_prev_even, vol_prev_odd, vol_merged
        type(image)                   :: vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd
        type(string)                  :: volname, eonames(2)
        type(restore_timings_t)       :: restore_timings
        type(nu_state_filter_timings) :: nu_timings
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: state_pops(:)
        logical, allocatable          :: l_state_dropped(:)
        real, allocatable             :: res0143s(:), res05s(:)
        real, allocatable             :: nu_align_lps(:)
        real, allocatable             :: update_frac_trail_recs(:), realized_update_fracs(:)
        integer                       :: state, numlen_part
        integer(timer_int_kind)       :: t_tot
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
            allocate(res0143s(params%nstates), res05s(params%nstates))
            res0143s = 0.
            res05s   = 0.
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
            res05s(state)           = 0.
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
                &volname, eonames, res05s(state), res0143s(state), restore_timings)
            params%vols(state)      = volname
            params%vols_even(state) = eonames(1)
            params%vols_odd(state)  = eonames(2)
            call postprocess_state()
            call volname%kill
        end subroutine assemble_state

        subroutine postprocess_state()
            if( l_nonuniform_mode ) call run_state_nonuniform_filter()
        end subroutine postprocess_state

        subroutine run_state_nonuniform_filter()
            ! the shared assembly-owned NU competition (both backends)
            call nonuniform_filter_state(params, state, vol_nu_base_even, vol_nu_base_odd, &
                &vol_nu_aux_even, vol_nu_aux_odd, use_static_nu_aux_replacement(), res0143s(state), &
                &volname, eonames, nu_align_lps(state), nu_timings)
        end subroutine run_state_nonuniform_filter

        logical function use_static_nu_aux_replacement() result(l_use_aux)
            l_use_aux = nu_static_aux_replacement(params)
        end function use_static_nu_aux_replacement

        subroutine collect_restore_timings()
            rt_nu_envmask                = nu_timings%envmask
            rt_nonuniform_filter         = nu_timings%filter
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
                call build%spproj_field%set_all2single('res',   res0143s(1))
                call build%spproj_field%set_all2single('res05', res05s(1))
            else
                do iptcl = 1, build%spproj_field%get_noris()
                    istate = build%spproj_field%get_state(iptcl)
                    if( istate > 0 .and. istate <= params%nstates )then
                        call build%spproj_field%set(iptcl, 'res',   res0143s(istate))
                        call build%spproj_field%set(iptcl, 'res05', res05s(istate))
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
            ! 'lp' is the raw handoff the matcher clips against lpstop;
            ! 'lp_est' records the NU-estimated limit itself for reporting
            call build%spproj_field%set_all2single('lp',     align_lp)
            call build%spproj_field%set_all2single('lp_est', align_lp)
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
            if( allocated(state_pops)             ) deallocate(state_pops)
            if( allocated(l_state_dropped)        ) deallocate(l_state_dropped)
            if( allocated(res0143s)               ) deallocate(res0143s)
            if( allocated(res05s)                 ) deallocate(res05s)
            if( allocated(nu_align_lps)           ) deallocate(nu_align_lps)
            if( allocated(update_frac_trail_recs) ) deallocate(update_frac_trail_recs)
            if( allocated(realized_update_fracs)  ) deallocate(realized_update_fracs)
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
