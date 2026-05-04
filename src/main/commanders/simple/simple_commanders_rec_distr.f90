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
    real(timer_int_kind) :: read_eos                       = 0.
    real(timer_int_kind) :: read_previous_halfmaps          = 0.
    real(timer_int_kind) :: sum_reduce                     = 0.
    real(timer_int_kind) :: sum_eos                        = 0.
    real(timer_int_kind) :: calc_fsc4sampl_dens_correct    = 0.
    real(timer_int_kind) :: sampl_dens_correct_eos         = 0.
    real(timer_int_kind) :: sampl_dens_correct_sum         = 0.
    real(timer_int_kind) :: insert_lowres                  = 0.
    real(timer_int_kind) :: trail_restored_halves          = 0.
end type restore_timings_t

contains

    !> Reduce one state's Cartesian partial reconstructions and restore dense
    !> even, odd, and merged volumes. On return build%vol/build%vol2 contain the
    !> restored half-volumes needed by automask and nonuniform filtering.
    subroutine restore_state_from_parts( params, build, cline, eorecvol_read, state, numlen_part, &
        &update_frac_trail_rec, gridcorr_img, vol_prev_even, vol_prev_odd, vol_e, vol_merged, &
        &volname, eonames, res05, res0143, timings )
        use simple_reconstructor_eo, only: reconstructor_eo
        type(parameters),       intent(in)    :: params
        type(builder),          intent(inout) :: build
        class(cmdline),         intent(in)    :: cline
        type(reconstructor_eo), intent(inout) :: eorecvol_read
        integer,                intent(in)    :: state, numlen_part
        real,                   intent(in)    :: update_frac_trail_rec
        type(image),            intent(inout) :: gridcorr_img
        type(image),            intent(inout) :: vol_prev_even, vol_prev_odd, vol_e, vol_merged
        type(string),           intent(inout) :: volname, eonames(2)
        real,                   intent(out)   :: res05, res0143
        type(restore_timings_t), intent(inout) :: timings
        type(string) :: volname_prev, volname_prev_even, volname_prev_odd, fsc_txt_file
        real, allocatable :: fsc(:)
        real    :: weight_prev
        integer :: find4eoavg
        integer(timer_int_kind) :: t_read_eos, t_read_previous_halfmaps, t_sum_reduce, t_sum_eos, t_calc_fsc
        integer(timer_int_kind) :: t_sampl_dens_correct, t_insert_lowres, t_trail
        call reduce_partials()
        call set_state_filenames()
        call sum_eos_before_density_correction_if_needed()
        call restore_eos_and_write_fsc()
        call sum_eos_after_density_correction_if_needed()
        call restore_merged_volume()
        ! Keep restored half-volumes current in build%vol/build%vol2 for
        ! automasking, nonuniform filtering, and optional trailing. In lp-set
        ! mode the halves are only on disk after density correction. Otherwise
        ! lowres_insert_into_halves leaves the odd half in build%vol2 and keeps
        ! a copy of the even half in vol_e while build%vol still holds the
        ! merged volume.
        if( params%l_lpset )then
            call build%vol%read(eonames(1))
            call build%vol2%read(eonames(2))
        else
            call build%vol%copy(vol_e)
        endif
        call trail_restored_halves_if_needed()
        call cleanup_restore_state()

    contains

        subroutine reduce_partials()
            integer :: part
            call build%eorecvol%reset_all
            do part = 1, params%nparts
                if( L_BENCH_GLOB ) t_read_eos = tic()
                call eorecvol_read%read_eos(refine3D_partial_rec_fbody(state, part, numlen_part))
                if( L_BENCH_GLOB )then
                    timings%read_eos = timings%read_eos + toc(t_read_eos)
                    t_sum_reduce = tic()
                endif
                call build%eorecvol%sum_reduce(eorecvol_read)
                if( L_BENCH_GLOB ) timings%sum_reduce = timings%sum_reduce + toc(t_sum_reduce)
            enddo
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
            if( params%l_trail_rec )then
                if( L_BENCH_GLOB ) t_read_previous_halfmaps = tic()
                call read_previous_halfmaps()
                if( L_BENCH_GLOB )then
                    timings%read_previous_halfmaps = timings%read_previous_halfmaps + toc(t_read_previous_halfmaps)
                endif
                if( allocated(fsc) ) deallocate(fsc)
                if( L_BENCH_GLOB ) t_calc_fsc = tic()
                call build%eorecvol%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                if( L_BENCH_GLOB )then
                    timings%calc_fsc4sampl_dens_correct = timings%calc_fsc4sampl_dens_correct + toc(t_calc_fsc)
                    t_sampl_dens_correct = tic()
                endif
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg, fsc)
            else
                if( L_BENCH_GLOB ) t_sampl_dens_correct = tic()
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            endif
            if( L_BENCH_GLOB )then
                timings%sampl_dens_correct_eos = timings%sampl_dens_correct_eos + toc(t_sampl_dens_correct)
            endif
            fsc_txt_file = resolve_fsc_txt_fname()
            call build%eorecvol%write_fsc2txt(fsc_txt_file)
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
            if( cline%defined('which_iter') )then
                fname = refine3D_resolution_txt_fbody(state, params%which_iter)
            else
                fname = refine3D_resolution_txt_fbody(state)
            endif
        end function resolve_fsc_txt_fname

        subroutine sum_eos_after_density_correction_if_needed()
            if( .not. params%l_ml_reg ) return
            if( L_BENCH_GLOB ) t_sum_eos = tic()
            call build%eorecvol%sum_eos
            if( L_BENCH_GLOB ) timings%sum_eos = timings%sum_eos + toc(t_sum_eos)
        end subroutine sum_eos_after_density_correction_if_needed

        subroutine restore_merged_volume()
            call build%eorecvol%get_res(res05, res0143)
            if( L_BENCH_GLOB ) t_sampl_dens_correct = tic()
            call build%eorecvol%sampl_dens_correct_sum(build%vol)
            if( L_BENCH_GLOB )then
                timings%sampl_dens_correct_sum = timings%sampl_dens_correct_sum + toc(t_sampl_dens_correct)
            endif
            call build%vol%fft
            if( .not. params%l_lpset )then
                if( L_BENCH_GLOB ) t_insert_lowres = tic()
                call lowres_insert_into_halves()
                if( L_BENCH_GLOB ) timings%insert_lowres = timings%insert_lowres + toc(t_insert_lowres)
            endif
            call build%vol%ifft
            call build%vol%mul(gridcorr_img)
            call build%vol%write(volname, del_if_exists=.true.)
            call wait_for_closure(volname)
        end subroutine restore_merged_volume

        subroutine lowres_insert_into_halves()
            call build%vol2%zero_and_unflag_ft
            call build%vol2%read(eonames(1))
            call build%vol2%fft()
            call build%vol2%insert_lowres(build%vol, find4eoavg)
            call build%vol2%ifft()
            call build%vol2%mul(gridcorr_img)
            call build%vol2%write(eonames(1), del_if_exists=.true.)
            call vol_e%copy(build%vol2)
            call build%vol2%zero_and_unflag_ft
            call build%vol2%read(eonames(2))
            call build%vol2%fft()
            call build%vol2%insert_lowres(build%vol, find4eoavg)
            call build%vol2%ifft()
            call build%vol2%mul(gridcorr_img)
            call build%vol2%write(eonames(2), del_if_exists=.true.)
        end subroutine lowres_insert_into_halves

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
            &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity
        use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess, AUTOMASK_ACTION_REGENERATE, &
            &NU_MASK_SOURCE_FRESH_AUTOMASK, NU_MASK_SOURCE_EXISTING_AUTOMASK
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(reconstructor_eo)        :: eorecvol_read
        type(image)                   :: vol_prev_even, vol_prev_odd, gridcorr_img, vol_merged
        type(image)                   :: vol_even_nu, vol_odd_nu, vol_msk
        type(image)                   :: vol_e, vol_nu_base_even, vol_nu_base_odd
        type(image), allocatable      :: nu_aux_even(:), nu_aux_odd(:)
        type(image_msk)               :: mskvol
        type(image_bin)               :: state_mask_bin
        type(string)                  :: volname, eonames(2)
        type(restore_timings_t)       :: restore_timings
        type(vol_pproc_plan)          :: pp_plan
        logical, allocatable          :: l_mask(:,:,:)
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: imat(:,:,:)
        real, allocatable             :: res0143s(:)
        real                          :: update_frac_trail_rec, res05
        integer                       :: state, ldim(3), ldim_pd(3), numlen_part
        integer(timer_int_kind)       :: t_automask3D, t_setup_nu_dmats, t_optimize_nu, t_nu_filter
        integer(timer_int_kind)       :: t_write_nu_outputs, t_tot
        real(timer_int_kind)          :: rt_read_eos, rt_sum_reduce, rt_sum_eos
        real(timer_int_kind)          :: rt_read_previous_halfmaps
        real(timer_int_kind)          :: rt_calc_fsc4sampl_dens_correct, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_insert_lowres, rt_trail_restored_halves
        real(timer_int_kind)          :: rt_automask3D, rt_setup_nu_dmats, rt_optimize_nu_cutoff_finds
        real(timer_int_kind)          :: rt_nu_filter_vols, rt_write_nonuniform_outputs, rt_tot
        call initialize_bench_timers()
        call initialize_context()
        call determine_trailing_update_fraction()
        call prepare_grid_correction()
        do state = 1, params%nstates
            call assemble_state()
        enddo
        call collect_restore_timings()
        call update_project_resolution_metadata()
        call cleanup_context()
        call write_benchmark()

    contains

        subroutine initialize_bench_timers()
            rt_read_eos                    = 0.
            rt_read_previous_halfmaps       = 0.
            rt_sum_reduce                  = 0.
            rt_sum_eos                     = 0.
            rt_calc_fsc4sampl_dens_correct = 0.
            rt_sampl_dens_correct_eos      = 0.
            rt_sampl_dens_correct_sum      = 0.
            rt_insert_lowres               = 0.
            rt_trail_restored_halves       = 0.
            rt_automask3D                  = 0.
            rt_setup_nu_dmats              = 0.
            rt_optimize_nu_cutoff_finds    = 0.
            rt_nu_filter_vols              = 0.
            rt_write_nonuniform_outputs     = 0.
            rt_tot                         = 0.
            if( L_BENCH_GLOB ) t_tot = tic()
        end subroutine initialize_bench_timers

        subroutine initialize_context()
            call build%init_params_and_build_general_tbox(cline, params)
            call build%build_rec_eo_tbox(params)
            call build%eorecvol%kill_exp
            numlen_part       = max(1, params%numlen)
            l_nonuniform_mode = trim(params%filt_mode).eq.'nonuniform'
            allocate(res0143s(params%nstates))
            res0143s = 0.
            call eorecvol_read%new(params, build%spproj, expand=.false.)
        end subroutine initialize_context

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
                &update_frac_trail_rec, gridcorr_img, vol_prev_even, vol_prev_odd, vol_e, vol_merged, &
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
            call build_nonuniform_mask()
            call setup_nonuniform_filter()
            if( L_BENCH_GLOB ) t_optimize_nu = tic()
            call optimize_nu_cutoff_finds()
            if( L_BENCH_GLOB ) rt_optimize_nu_cutoff_finds = rt_optimize_nu_cutoff_finds + toc(t_optimize_nu)
            if( L_BENCH_GLOB ) t_nu_filter = tic()
            call nu_filter_vols(vol_even_nu, vol_odd_nu)
            if( L_BENCH_GLOB ) rt_nu_filter_vols = rt_nu_filter_vols + toc(t_nu_filter)
            call log_nonuniform_filter_stats()
            if( L_BENCH_GLOB ) t_write_nu_outputs = tic()
            call write_nonuniform_outputs()
            if( L_BENCH_GLOB )then
                rt_write_nonuniform_outputs = rt_write_nonuniform_outputs + toc(t_write_nu_outputs)
            endif
            call cleanup_nonuniform_state()
        end subroutine run_state_nonuniform_filter

        subroutine build_nonuniform_mask()
            real :: mskrad_px
            if( allocated(l_mask) ) deallocate(l_mask)
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
            if( allocated(nu_aux_even) ) deallocate(nu_aux_even)
            if( allocated(nu_aux_odd) )  deallocate(nu_aux_odd)
            if( params%l_ml_reg )then
                call vol_nu_base_even%new(ldim, params%smpd_crop)
                call vol_nu_base_odd%new( ldim, params%smpd_crop)
                call vol_nu_base_even%read(add2fbody(eonames(1), MRC_EXT, '_unfil'))
                call vol_nu_base_odd%read( add2fbody(eonames(2), MRC_EXT, '_unfil'))
                allocate(nu_aux_even(1), nu_aux_odd(1))
                call nu_aux_even(1)%copy(build%vol)
                call nu_aux_odd(1)%copy(build%vol2)
                if( L_BENCH_GLOB ) t_setup_nu_dmats = tic()
                call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, l_mask, [res0143s(state)], &
                    &nu_aux_even, nu_aux_odd)
            else
                if( L_BENCH_GLOB ) t_setup_nu_dmats = tic()
                call setup_nu_dmats(build%vol, build%vol2, l_mask, [real ::])
            endif
            if( L_BENCH_GLOB ) rt_setup_nu_dmats = rt_setup_nu_dmats + toc(t_setup_nu_dmats)
        end subroutine setup_nonuniform_filter

        subroutine log_nonuniform_filter_stats()
            if( allocated(nu_aux_even) )then
                call print_nu_filtmap_lowpass_stats(l_mask, aux_resolutions=[res0143s(state)])
            else
                call print_nu_filtmap_lowpass_stats(l_mask)
            endif
            call analyze_filtmap_neighbor_continuity(l_mask)
        end subroutine log_nonuniform_filter_stats

        subroutine write_nonuniform_outputs()
            type(string) :: eonames_nu(2), volname_nu
            eonames_nu(1) = add2fbody(eonames(1), MRC_EXT, NUFILT_SUFFIX)
            eonames_nu(2) = add2fbody(eonames(2), MRC_EXT, NUFILT_SUFFIX)
            volname_nu    = add2fbody(volname,    MRC_EXT, NUFILT_SUFFIX)
            call vol_even_nu%write(eonames_nu(1), del_if_exists=.true.)
            call vol_odd_nu%write(eonames_nu(2), del_if_exists=.true.)
            call vol_even_nu%add(vol_odd_nu)
            call vol_even_nu%mul(0.5)
            call vol_even_nu%write(volname_nu, del_if_exists=.true.)
            call wait_for_closure(volname_nu)
            call eonames_nu(1)%kill
            call eonames_nu(2)%kill
            call volname_nu%kill
        end subroutine write_nonuniform_outputs

        subroutine cleanup_nonuniform_state()
            call vol_even_nu%kill
            call vol_odd_nu%kill
            call vol_nu_base_even%kill
            call vol_nu_base_odd%kill
            if( allocated(nu_aux_even) )then
                call nu_aux_even(1)%kill
                deallocate(nu_aux_even)
            endif
            if( allocated(nu_aux_odd) )then
                call nu_aux_odd(1)%kill
                deallocate(nu_aux_odd)
            endif
            call vol_msk%kill
            if( allocated(l_mask) ) deallocate(l_mask)
            call cleanup_nu_filter()
        end subroutine cleanup_nonuniform_state

        subroutine collect_restore_timings()
            rt_read_eos                    = restore_timings%read_eos
            rt_read_previous_halfmaps       = restore_timings%read_previous_halfmaps
            rt_sum_reduce                  = restore_timings%sum_reduce
            rt_sum_eos                     = restore_timings%sum_eos
            rt_calc_fsc4sampl_dens_correct = restore_timings%calc_fsc4sampl_dens_correct
            rt_sampl_dens_correct_eos      = restore_timings%sampl_dens_correct_eos
            rt_sampl_dens_correct_sum      = restore_timings%sampl_dens_correct_sum
            rt_insert_lowres               = restore_timings%insert_lowres
            rt_trail_restored_halves       = restore_timings%trail_restored_halves
        end subroutine collect_restore_timings

        subroutine update_project_resolution_metadata()
            integer :: iptcl, istate
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
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        end subroutine update_project_resolution_metadata

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
            if( allocated(nu_aux_even) )then
                call nu_aux_even(1)%kill
                deallocate(nu_aux_even)
            endif
            if( allocated(nu_aux_odd) )then
                call nu_aux_odd(1)%kill
                deallocate(nu_aux_odd)
            endif
            if( allocated(l_mask) ) deallocate(l_mask)
            if( allocated(imat) ) deallocate(imat)
            call cleanup_nu_filter()
            call state_mask_bin%kill_bimg
            call mskvol%kill_bimg
            call vol_e%kill
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
            write(fnr,'(a)')    'volassemble assembly mode           : volume'
            write(fnr,'(a,i0)') 'volassemble nspace                  : ', params%nspace
            write(fnr,'(a,i0)') 'volassemble nstates                 : ', params%nstates
            write(fnr,'(a,i0)') 'volassemble kfrom                   : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'volassemble kto                     : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble read_eos                : ', rt_read_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble read_previous_halfmaps  : ', rt_read_previous_halfmaps
            write(fnr,'(a,t52,f9.2)') 'volassemble sum_reduce              : ', rt_sum_reduce
            write(fnr,'(a,t52,f9.2)') 'volassemble sum_eos                 : ', rt_sum_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble calc_fsc4sampl_dens_correct : ', &
                rt_calc_fsc4sampl_dens_correct
            write(fnr,'(a,t52,f9.2)') 'volassemble sampl_dens_correct_eos  : ', rt_sampl_dens_correct_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble sampl_dens_correct_sum  : ', rt_sampl_dens_correct_sum
            write(fnr,'(a,t52,f9.2)') 'volassemble insert_lowres           : ', rt_insert_lowres
            write(fnr,'(a,t52,f9.2)') 'volassemble trail_restored_halves   : ', rt_trail_restored_halves
            write(fnr,'(a,t52,f9.2)') 'volassemble automask3D              : ', rt_automask3D
            write(fnr,'(a,t52,f9.2)') 'volassemble setup_nu_dmats          : ', rt_setup_nu_dmats
            write(fnr,'(a,t52,f9.2)') 'volassemble optimize_nu_cutoff_finds : ', rt_optimize_nu_cutoff_finds
            write(fnr,'(a,t52,f9.2)') 'volassemble nu_filter_vols          : ', rt_nu_filter_vols
            write(fnr,'(a,t52,f9.2)') 'volassemble write_nonuniform_outputs : ', rt_write_nonuniform_outputs
            write(fnr,'(a,t52,f9.2)') 'volassemble total time              : ', rt_tot
            call fclose(fnr)
            call benchfname%kill
        end subroutine write_benchmark

    end subroutine exec_volassemble

end module simple_commanders_rec_distr
