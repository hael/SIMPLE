module simple_commanders_rec_distr
use simple_commanders_api
use simple_refine3D_fnames, only: refine3D_partial_rec_fbody, refine3D_resolution_txt_fbody, &
    &refine3D_state_halfvol_fname, refine3D_state_vol_fname, refine3D_volassemble_bench_fname
implicit none
private
public :: commander_volassemble
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
    real(timer_int_kind) :: trail_restored_halves          = 0.
end type restore_timings_t

contains

    !> Reduce one state's Cartesian partial reconstructions and restore dense
    !> even, odd, and merged volumes. On return build%vol/build%vol2 contain the
    !> restored half-volumes needed by automask, while vol_nu_base_even/odd and
    !> optional vol_nu_aux_even/odd contain the static-bank nonuniform-filter
    !> auxiliary inputs before even/odd low-resolution insertion.
    subroutine restore_state_from_parts( params, build, cline, eorecvol_read, state, numlen_part, &
        &update_frac_trail_rec, gridcorr_img, vol_prev_even, vol_prev_odd, vol_merged, &
        &vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd, &
        &volname, eonames, res05, res0143, timings )
        use simple_reconstructor_eo, only: reconstructor_eo
        type(parameters),       intent(in)    :: params
        type(builder),          intent(inout) :: build
        class(cmdline),         intent(in)    :: cline
        type(reconstructor_eo), intent(inout) :: eorecvol_read
        integer,                intent(in)    :: state, numlen_part
        real,                   intent(in)    :: update_frac_trail_rec
        type(image),            intent(inout) :: gridcorr_img
        type(image),            intent(inout) :: vol_prev_even, vol_prev_odd, vol_merged
        type(image),            intent(inout) :: vol_nu_base_even, vol_nu_base_odd
        type(image),            intent(inout) :: vol_nu_aux_even, vol_nu_aux_odd
        type(string),           intent(inout) :: volname, eonames(2)
        real,                   intent(out)   :: res05, res0143
        type(restore_timings_t), intent(inout) :: timings
        type(string) :: volname_prev, volname_prev_even, volname_prev_odd
        type(string) :: fsc_txt_file
        real, allocatable :: fsc(:)
        real    :: weight_prev
        integer :: find4eoavg, ldim(3)
        integer(timer_int_kind) :: t_reduce_partials, t_restore_eos, t_restore_merged, t_sum_eos, t_trail
        call reduce_partials()
        call set_state_filenames()
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
            call build%eorecvol%reset_all
            do part = 1, params%nparts
                call eorecvol_read%read_eos(refine3D_partial_rec_fbody(state, part, numlen_part))
                call build%eorecvol%sum_reduce(eorecvol_read)
            enddo
            if( L_BENCH_GLOB ) timings%reduce_partials = timings%reduce_partials + toc(t_reduce_partials)
        end subroutine reduce_partials

        subroutine set_state_filenames()
            volname    = refine3D_state_vol_fname(state)
            eonames(1) = refine3D_state_halfvol_fname(state, 'even')
            eonames(2) = refine3D_state_halfvol_fname(state, 'odd')
        end subroutine set_state_filenames

        subroutine sum_eos_before_density_correction_if_needed()
            if( params%l_ml_reg ) return
            if( L_BENCH_GLOB ) t_sum_eos = tic()
            call build%eorecvol%sum_eos
            if( L_BENCH_GLOB ) timings%sum_eos = timings%sum_eos + toc(t_sum_eos)
        end subroutine sum_eos_before_density_correction_if_needed

        subroutine restore_eos_and_write_fsc()
            if( L_BENCH_GLOB ) t_restore_eos = tic()
            if( params%l_trail_rec )then
                call read_previous_halfmaps()
                if( allocated(fsc) ) deallocate(fsc)
                call build%eorecvol%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg, fsc)
            else
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            endif
            fsc_txt_file = resolve_fsc_txt_fname()
            call build%eorecvol%write_fsc2txt(fsc_txt_file)
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
            if( L_BENCH_GLOB ) t_sum_eos = tic()
            call build%eorecvol%sum_eos
            if( L_BENCH_GLOB ) timings%sum_eos = timings%sum_eos + toc(t_sum_eos)
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
                    call vol_nu_aux_even%mul(gridcorr_img)
                    call vol_nu_aux_odd%mul( gridcorr_img)
                endif
            else
                call vol_nu_base_even%read(eonames(1))
                call vol_nu_base_odd%read( eonames(2))
            endif
            call vol_nu_base_even%mul(gridcorr_img)
            call vol_nu_base_odd%mul( gridcorr_img)
        end subroutine capture_nonuniform_source_halves

        subroutine restore_merged_volume()
            if( L_BENCH_GLOB ) t_restore_merged = tic()
            call build%eorecvol%get_res(res05, res0143)
            call build%eorecvol%sampl_dens_correct_sum(build%vol)
            call build%vol%fft
            call build%vol%ifft
            call build%vol%mul(gridcorr_img)
            call build%vol%write(volname, del_if_exists=.true.)
            call wait_for_closure(volname)
            if( L_BENCH_GLOB ) timings%restore_merged_volume = &
                timings%restore_merged_volume + toc(t_restore_merged)
        end subroutine restore_merged_volume

        subroutine trail_restored_halves_if_needed()
            if( .not. params%l_trail_rec ) return
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
            call volname_prev%kill
            call volname_prev_even%kill
            call volname_prev_odd%kill
            if( allocated(fsc) ) deallocate(fsc)
        end subroutine cleanup_restore_state

    end subroutine restore_state_from_parts

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor_eo, only: reconstructor_eo
        use simple_gridding,         only: prep3D_inv_instrfun4mul
        use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
            &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity, &
            &extend_nu_filter_highres_shell_next, refine_nu_extension_filtmap_ordered_labels, &
            &nu_highres_extension_stats, get_nu_filtmap_finest_selected_lp, &
            &get_nu_filtmap_highres_shell_depth, write_nu_local_resolution_map
        use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess, AUTOMASK_ACTION_REGENERATE, &
            &NU_MASK_SOURCE_FRESH_AUTOMASK, NU_MASK_SOURCE_EXISTING_AUTOMASK
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(reconstructor_eo)        :: eorecvol_read
        type(image)                   :: vol_prev_even, vol_prev_odd, gridcorr_img, vol_merged
        type(image)                   :: vol_even_nu, vol_odd_nu, vol_msk
        type(image)                   :: vol_nu_base_even, vol_nu_base_odd, vol_nu_aux_even, vol_nu_aux_odd
        type(image), allocatable      :: nu_aux_even(:), nu_aux_odd(:)
        type(image_msk)               :: mskvol
        type(image_bin)               :: state_mask_bin
        type(string)                  :: volname, eonames(2)
        type(restore_timings_t)       :: restore_timings
        type(vol_pproc_plan)          :: pp_plan
        logical, allocatable          :: l_mask(:,:,:)
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: imat(:,:,:)
        integer, allocatable          :: state_pops(:)
        real, allocatable             :: res0143s(:)
        real, allocatable             :: nu_align_lps(:)
        real                          :: update_frac_trail_rec, res05
        integer                       :: state, ldim(3), ldim_pd(3), numlen_part
        integer(timer_int_kind)       :: t_automask3D, t_nonuniform_filter, t_tot
        integer(timer_int_kind)       :: t_init_context, t_trail_frac, t_gridcorr, t_upd_proj, t_cleanup
        real(timer_int_kind)          :: rt_reduce_partials, rt_sum_eos
        real(timer_int_kind)          :: rt_restore_eos_and_write_fsc, rt_restore_merged_volume
        real(timer_int_kind)          :: rt_trail_restored_halves
        real(timer_int_kind)          :: rt_automask3D, rt_nonuniform_filter, rt_tot
        real(timer_int_kind)          :: rt_init_context, rt_trail_frac, rt_gridcorr, rt_upd_proj, rt_cleanup
        call initialize_bench_timers()
        if( L_BENCH_GLOB ) t_init_context = tic()
        call initialize_context()
        if( L_BENCH_GLOB ) rt_init_context = toc(t_init_context)
        if( L_BENCH_GLOB ) t_trail_frac = tic()
        call determine_trailing_update_fraction()
        if( L_BENCH_GLOB ) rt_trail_frac = toc(t_trail_frac)
        if( L_BENCH_GLOB ) t_gridcorr = tic()
        call prepare_grid_correction()
        if( L_BENCH_GLOB ) rt_gridcorr = toc(t_gridcorr)
        do state = 1, params%nstates
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
            rt_trail_restored_halves     = 0.
            rt_automask3D                = 0.
            rt_nonuniform_filter       = 0.
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
            call build%build_rec_eo_tbox(params)
            call build%eorecvol%kill_exp
            numlen_part       = max(1, params%numlen)
            l_nonuniform_mode = params%l_nonuniform
            allocate(res0143s(params%nstates))
            res0143s = 0.
            allocate(nu_align_lps(params%nstates))
            nu_align_lps = 0.
            allocate(state_pops(params%nstates))
            state_pops = 0
            call eorecvol_read%new(params, build%spproj, expand=.false.)
        end subroutine initialize_context

        subroutine refresh_state_populations()
            integer :: istate
            if( .not. allocated(state_pops) ) return
            do istate = 1,params%nstates
                state_pops(istate) = build%spproj_field%get_pop(istate, 'state')
            enddo
        end subroutine refresh_state_populations

        subroutine determine_trailing_update_fraction()
            update_frac_trail_rec = 1.0
            if( .not. params%l_trail_rec ) return
            if( cline%defined('ufrac_trec') )then
                update_frac_trail_rec = params%ufrac_trec
            else
                call build%spproj%read_segment(params%oritype, params%projfile)
                update_frac_trail_rec = build%spproj%os_ptcl3D%get_update_frac()
            endif
        end subroutine determine_trailing_update_fraction

        subroutine prepare_grid_correction()
            ldim         = build%vol%get_ldim()
            ldim_pd      = OSMPL_PAD_FAC * ldim
            gridcorr_img = prep3D_inv_instrfun4mul(ldim, ldim_pd, params%smpd_crop)
        end subroutine prepare_grid_correction

        subroutine assemble_state()
            call restore_state_from_parts(params, build, cline, eorecvol_read, state, numlen_part, &
                &update_frac_trail_rec, gridcorr_img, vol_prev_even, vol_prev_odd, vol_merged, &
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
            call plan_state_postprocess(params, state, which_iter, l_nonuniform_mode, pp_plan)
            if( pp_plan%l_state_mask_incompatible )then
                write(logfhandle,'(A,1X,A)') &
                    &'>>> Existing automask incompatible with current box/sampling, regenerating:', &
                    &pp_plan%mskfile_state%to_char()
            endif
            call run_state_automask()
            if( l_nonuniform_mode ) call run_state_nonuniform_filter()
        end subroutine postprocess_state

        subroutine run_state_automask()
            if( pp_plan%automask_action == AUTOMASK_ACTION_REGENERATE )then
                if( L_BENCH_GLOB ) t_automask3D = tic()
                call mskvol%automask3D(params, build%vol, build%vol2, pp_plan%automask_tight)
                if( L_BENCH_GLOB ) rt_automask3D = rt_automask3D + toc(t_automask3D)
                call mskvol%write(pp_plan%mskfile_state)
            endif
        end subroutine run_state_automask

        subroutine run_state_nonuniform_filter()
            if( L_BENCH_GLOB ) t_nonuniform_filter = tic()
            call build_nonuniform_mask()
            call setup_nonuniform_filter()
            if( allocated(l_mask) ) deallocate(l_mask)
            call release_nonuniform_aux_inputs()
            call optimize_nu_cutoff_finds()
            call refine_nonuniform_filter_bank()
            call release_nonuniform_base_inputs()
            call nu_filter_vols(vol_even_nu, vol_odd_nu, soft_synthesis=params%l_nu_soft_synth)
            call log_nonuniform_filter_stats()
            call write_nonuniform_outputs()
            call record_nu_alignment_lowpass_limit()
            call cleanup_nonuniform_state()
            if( L_BENCH_GLOB ) rt_nonuniform_filter = rt_nonuniform_filter + toc(t_nonuniform_filter)
        end subroutine run_state_nonuniform_filter

        subroutine build_nonuniform_mask()
            real :: mskrad_px
            if( allocated(l_mask) ) deallocate(l_mask)
            if( allocated(imat) ) deallocate(imat)
            select case( pp_plan%nu_mask_source )
                case( NU_MASK_SOURCE_FRESH_AUTOMASK )
                    call mskvol%set_imat
                    call mskvol%get_imat(imat)
                case( NU_MASK_SOURCE_EXISTING_AUTOMASK )
                    call state_mask_bin%new_bimg(ldim, params%smpd_crop)
                    call state_mask_bin%read_bimg(pp_plan%mskfile_state)
                    call state_mask_bin%get_imat(imat)
                    call state_mask_bin%kill_bimg
            end select
            if( allocated(imat) )then
                allocate(l_mask(ldim(1),ldim(2),ldim(3)))
                l_mask = imat > 0
                deallocate(imat)
            endif
            if( .not. allocated(l_mask) )then
                mskrad_px = 0.5 * params%mskdiam / params%smpd_crop
                call vol_msk%disc(ldim, params%smpd_crop, mskrad_px, l_mask)
            endif
        end subroutine build_nonuniform_mask

        subroutine setup_nonuniform_filter()
            integer :: n_highres_steps
            n_highres_steps = nu_highres_steps_for_state()
            call cleanup_nu_aux_images()
            if( use_static_nu_aux_replacement() )then
                allocate(nu_aux_even(1), nu_aux_odd(1))
                call nu_aux_even(1)%copy(vol_nu_aux_even)
                call nu_aux_odd(1)%copy(vol_nu_aux_odd)
                call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, l_mask, [res0143s(state)], &
                    &nu_aux_even, nu_aux_odd, n_highres_steps=n_highres_steps)
            else
                call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, l_mask, [real ::], &
                    &n_highres_steps=n_highres_steps)
            endif
        end subroutine setup_nonuniform_filter

        logical function use_static_nu_aux_replacement() result(l_use_aux)
            l_use_aux = params%l_ml_reg .and. .not. params%l_nu_refine
        end function use_static_nu_aux_replacement

        subroutine refine_nonuniform_filter_bank()
            type(nu_highres_extension_stats) :: ext_stats
            integer :: n_highres_steps, n_accepted_this_iteration
            if( .not. params%l_nu_refine ) return
            n_highres_steps = nu_highres_steps_for_state()
            n_accepted_this_iteration = 0
            do
                call extend_nu_filter_highres_shell_next(vol_nu_base_even, vol_nu_base_odd, stats=ext_stats)
                if( .not. ext_stats%attempted )then
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
                    exit
                endif
                if( .not. ext_stats%applied      ) exit
                if( .not. ext_stats%promote_next ) exit
                n_accepted_this_iteration = n_accepted_this_iteration + 1
            end do
            if( n_accepted_this_iteration > 0 )then
                call refine_nu_extension_filtmap_ordered_labels()
                n_highres_steps = get_nu_filtmap_highres_shell_depth()
                call write_nu_highres_steps_for_state(n_highres_steps)
                write(logfhandle,'(A,I0,A,I0)') &
                    &'>>> NU high-resolution extension accepted shell steps this iteration: ', &
                    &n_accepted_this_iteration, '; promoted depth for next iteration: ', n_highres_steps
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
            call analyze_filtmap_neighbor_continuity()
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
            real :: align_lp
            if( .not. params%l_nonuniform ) return
            align_lp = get_nu_filtmap_finest_selected_lp()
            if( align_lp <= TINY ) return
            nu_align_lps(state) = align_lp
            write(logfhandle,'(A,I0,A,F8.3,A)') &
                &'>>> NU filter state ', state, ' matching low-pass limit for next iteration: ', align_lp, ' A'
        end subroutine record_nu_alignment_lowpass_limit

        subroutine cleanup_nonuniform_state()
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
            call vol_nu_aux_even%kill
            call vol_nu_aux_odd%kill
            call cleanup_nu_aux_images()
            call vol_msk%kill
            if( allocated(l_mask) ) deallocate(l_mask)
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
            if( params%nstates > 1 ) call log_nu_alignment_lowpass_summary(l_included)
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
            write(logfhandle,'(A,I0,A,F8.3,A)') &
                &'>>> NU filter project matching low-pass limit from state ', selected_state, ': ', align_lp, ' A'
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
            call gridcorr_img%kill
            call build%kill_general_tbox
            call build%eorecvol%kill_exp
            call build%kill_rec_eo_tbox
            call eorecvol_read%kill
            call vol_prev_even%kill
            call vol_prev_odd%kill
            call vol_merged%kill
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
            call vol_nu_aux_even%kill
            call vol_nu_aux_odd%kill
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call vol_msk%kill
            call cleanup_nu_aux_images()
            if( allocated(l_mask) ) deallocate(l_mask)
            if( allocated(imat) ) deallocate(imat)
            if( allocated(state_pops) ) deallocate(state_pops)
            if( allocated(res0143s) ) deallocate(res0143s)
            if( allocated(nu_align_lps) ) deallocate(nu_align_lps)
            call cleanup_nu_filter()
            call state_mask_bin%kill_bimg
            call mskvol%kill_bimg
            call pp_plan%mskfile_state%kill
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
            write(fnr,'(a,t52,f9.2)') 'volassemble reduce_partials           : ', rt_reduce_partials
            write(fnr,'(a,t52,f9.2)') 'volassemble sum_eos                   : ', rt_sum_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble restore_eos_and_write_fsc : ', &
                rt_restore_eos_and_write_fsc
            write(fnr,'(a,t52,f9.2)') 'volassemble restore_merged_volume     : ', rt_restore_merged_volume
            write(fnr,'(a,t52,f9.2)') 'volassemble trail_restored_halves     : ', rt_trail_restored_halves
            write(fnr,'(a,t52,f9.2)') 'volassemble automask3D                : ', rt_automask3D
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform_filter         : ', rt_nonuniform_filter
            write(fnr,'(a,t52,f9.2)') 'volassemble init_context              : ', rt_init_context
            write(fnr,'(a,t52,f9.2)') 'volassemble trail_frac_read           : ', rt_trail_frac
            write(fnr,'(a,t52,f9.2)') 'volassemble gridcorr_prep             : ', rt_gridcorr
            write(fnr,'(a,t52,f9.2)') 'volassemble update_proj_metadata      : ', rt_upd_proj
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup                   : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble total time                : ', rt_tot
            write(fnr,'(a,t52,f9.2)') 'volassemble % accounted for           : ', &
                &((rt_reduce_partials + rt_sum_eos + rt_restore_eos_and_write_fsc +               &
                &  rt_restore_merged_volume + rt_trail_restored_halves + rt_automask3D +          &
                &  rt_nonuniform_filter +                                                         &
                &  rt_init_context + rt_trail_frac + rt_gridcorr + rt_upd_proj + rt_cleanup)      &
                & / rt_tot) * 100.
            call fclose(fnr)
            call benchfname%kill
        end subroutine write_benchmark

    end subroutine exec_volassemble

end module simple_commanders_rec_distr
