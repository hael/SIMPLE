module simple_commanders_rec_distr
use simple_commanders_api
use simple_refine3D_fnames, only: refine3D_fsc_fname, refine3D_iter_refs_fname, refine3D_polar_refs_fbody, &
    &refine3D_partial_rec_fbody, refine3D_resolution_txt_fbody, refine3D_state_halfvol_fname, &
    &refine3D_state_vol_fbody, refine3D_state_vol_fname, refine3D_volassemble_bench_fname
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_cartesian_volassemble
  ! Cartesian assembly path.
  contains
    procedure :: execute => exec_cartesian_assembly
end type commander_cartesian_volassemble

type, extends(commander_base) :: commander_polar_volassemble
  ! Polar reference assembly path for polar=yes|obsfield.
  contains
    procedure :: execute => exec_polar_assembly
end type commander_polar_volassemble

type, extends(commander_base) :: commander_polar_rec3D_worker
  ! Polar partial-reference worker for distributed reconstruct3D.
  contains
    procedure :: execute => exec_polar_rec3D_worker
end type commander_polar_rec3D_worker

contains

    subroutine exec_polar_rec3D_worker( self, cline )
        use simple_matcher_3Dpolar,          only: calc_polar_partials
        use simple_matcher_smpl_and_lplims,  only: set_bp_range3D
        use simple_qsys_funs,                only: qsys_job_finished
        class(commander_polar_rec3D_worker), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(string)         :: fname
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( .not. params%l_polar ) THROW_HARD('polar_rec3D worker requires POLAR=yes|obsfield')
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        if( params%l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            call build%esig%read_groups(build%spproj_field)
        end if
        call set_bp_range3D(params, build, cline)
        call calc_polar_partials(params, build, cline, nptcls2update, pinds)
        if( allocated(pinds) ) deallocate(pinds)
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_rec_distr :: exec_polar_rec3D_worker'))
    end subroutine exec_polar_rec3D_worker

    subroutine exec_polar_assembly( self, cline )
        use simple_gridding,                only: prep3D_inv_instrfun4mul
        use simple_matcher_refvol_utils,    only: write_polar_refs_from_current_pftc
        use simple_matcher_smpl_and_lplims, only: set_bp_range3D
        use simple_polarft_calc,            only: vol_pad2ref_pfts
        use simple_reconstructor_eo,        only: reconstructor_eo
        use simple_class_frcs,              only: class_frcs
        class(commander_polar_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(reconstructor_eo) :: eorecvol_read
        type(class_frcs)       :: obs_frcs
        type(image)            :: vol_prev_even, vol_prev_odd, gridcorr_img, vol_e, vol_merged
        type(string)           :: recname, volname, volname_prev, fsc_file, fsc_txt_file
        type(string)           :: volname_prev_even, volname_prev_odd, eonames(2), benchfname
        real, allocatable      :: fsc(:), fsc_state(:)
        integer(timer_int_kind) :: t_tot, t_setup, t_reduce, t_normalize, t_resolution, t_write, t_cleanup
        integer(timer_int_kind) :: t_obs_detail
        real(timer_int_kind)    :: rt_tot, rt_setup, rt_reduce, rt_normalize, rt_resolution, rt_write, rt_cleanup
        real(timer_int_kind)    :: rt_cmp_restore, rt_cmp_accounted
        real(timer_int_kind)    :: rt_obs_reduce_detail(2), rt_obs_norm_detail(20)
        real                   :: update_frac_eff, update_frac_trail_rec, weight_prev
        integer :: fnr, nrefs, state, part, numlen_part, find4eoavg, ldim(3), ldim_pd(3), iproj
        logical :: l_refs_written
        rt_tot               = 0.
        rt_setup             = 0.
        rt_reduce            = 0.
        rt_normalize         = 0.
        rt_resolution        = 0.
        rt_write             = 0.
        rt_cleanup           = 0.
        rt_obs_reduce_detail = 0.
        rt_obs_norm_detail   = 0.
        l_refs_written       = .false.
        call build%init_params_and_build_general_tbox(cline, params)
        if( L_BENCH_GLOB ) t_tot = tic()
        ! Matchers write partition-local Cartesian partial reconstructions,
        ! polar partial sums (polar=yes), or obsfield Cartesian partials, while
        ! the assembly commander owns shared-memory reduction and reference update.
        if( L_BENCH_GLOB ) t_setup = tic()
        call set_bp_range3D(params, build, cline)
        ! The assembly cline may have promoted nspace_next/pftsz_next into the
        ! live nspace/pftsz values; project obsfield volumes on that grid.
        if( build%eulspace%get_noris() /= params%nspace )then
            call build%eulspace%kill
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
        endif
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,1], params%kfromto)
        params%refs = refine3D_iter_refs_fname(params%which_iter)
        call build%pftc%polar_cavger_new(.true., nrefs=nrefs)
        call build%pftc%polar_cavger_calc_pops(build%spproj)
        update_frac_eff = params%update_frac
        if( params%l_trail_rec .and. (.not. cline%defined('ufrac_trec')) )then
            if( build%spproj_field%has_been_sampled() )then
                update_frac_eff = build%spproj_field%get_update_frac()
            endif
        endif
        if( L_BENCH_GLOB ) rt_setup = toc(t_setup)
        select case(trim(params%polar))
            case('obsfield')
                call build%build_rec_eo_tbox(params)
                call eorecvol_read%new(params, build%spproj, expand=.false.)
                numlen_part = max(1, params%numlen)
                update_frac_trail_rec = 1.0
                if( params%l_trail_rec )then
                    if( cline%defined('ufrac_trec') )then
                        update_frac_trail_rec = params%ufrac_trec
                    else
                        update_frac_trail_rec = update_frac_eff
                    endif
                endif
                ldim         = build%vol%get_ldim()
                ldim_pd      = OSMPL_PAD_FAC * ldim
                gridcorr_img = prep3D_inv_instrfun4mul(ldim, ldim_pd, params%smpd_crop)
                call obs_frcs%new(params%nspace, params%box_crop, params%smpd_crop, params%nstates)
                do state = 1, params%nstates
                    call build%eorecvol%reset_all
                    do part = 1, params%nparts
                        if( L_BENCH_GLOB ) t_obs_detail = tic()
                        call eorecvol_read%read_eos(refine3D_partial_rec_fbody(state, part, numlen_part))
                        if( L_BENCH_GLOB )then
                            rt_obs_reduce_detail(1) = rt_obs_reduce_detail(1) + toc(t_obs_detail)
                            t_obs_detail = tic()
                        endif
                        call build%eorecvol%sum_reduce(eorecvol_read)
                        if( L_BENCH_GLOB ) rt_obs_reduce_detail(2) = rt_obs_reduce_detail(2) + toc(t_obs_detail)
                    enddo
                    if( L_BENCH_GLOB ) t_normalize = tic()
                    recname    = refine3D_state_vol_fbody(state)
                    volname    = refine3D_state_vol_fname(state)
                    eonames(1) = refine3D_state_halfvol_fname(state, 'even')
                    eonames(2) = refine3D_state_halfvol_fname(state, 'odd')
                    if( params%l_ml_reg )then
                        ! The Cartesian path sums even/odd after ML regularization.
                    else
                        call build%eorecvol%sum_eos
                    endif
                    if( params%l_trail_rec )then
                        if( cline%defined('vol'//int2str(state)) )then
                            volname_prev = cline%get_carg('vol'//int2str(state))
                        else
                            volname_prev = refine3D_state_vol_fname(state)
                        endif
                        volname_prev_even = add2fbody(volname_prev, MRC_EXT, '_even')
                        volname_prev_odd  = add2fbody(volname_prev, MRC_EXT, '_odd')
                        if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
                        if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
                        if( L_BENCH_GLOB ) t_obs_detail = tic()
                        call vol_prev_even%read_and_crop(volname_prev_even, params%smpd, params%box_crop, params%smpd_crop)
                        call vol_prev_odd%read_and_crop( volname_prev_odd,  params%smpd, params%box_crop, params%smpd_crop)
                        if( L_BENCH_GLOB ) rt_obs_norm_detail(4) = rt_obs_norm_detail(4) + toc(t_obs_detail)
                        if( allocated(fsc) ) deallocate(fsc)
                        if( L_BENCH_GLOB ) t_obs_detail = tic()
                        call build%eorecvol%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                        call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg, fsc)
                        if( L_BENCH_GLOB ) rt_obs_norm_detail(12) = rt_obs_norm_detail(12) + toc(t_obs_detail)
                    else
                        if( L_BENCH_GLOB ) t_obs_detail = tic()
                        call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
                        if( L_BENCH_GLOB ) rt_obs_norm_detail(12) = rt_obs_norm_detail(12) + toc(t_obs_detail)
                    endif
                    if( cline%defined('which_iter') )then
                        fsc_txt_file = refine3D_resolution_txt_fbody(state, params%which_iter)
                    else
                        fsc_txt_file = refine3D_resolution_txt_fbody(state)
                    endif
                    call build%eorecvol%write_fsc2txt(fsc_txt_file)
                    if( params%l_ml_reg ) call build%eorecvol%sum_eos
                    if( L_BENCH_GLOB ) t_obs_detail = tic()
                    call build%eorecvol%sampl_dens_correct_sum(build%vol)
                    if( L_BENCH_GLOB ) rt_obs_norm_detail(7) = rt_obs_norm_detail(7) + toc(t_obs_detail)
                    if( L_BENCH_GLOB ) t_obs_detail = tic()
                    call build%vol%fft
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
                    call build%vol%copy(vol_e)
                    call vol_merged%copy(vol_e)
                    call vol_merged%add(build%vol2)
                    call vol_merged%mul(0.5)
                    call vol_merged%write(volname, del_if_exists=.true.)
                    call wait_for_closure(volname)
                    if( L_BENCH_GLOB ) rt_obs_norm_detail(6) = rt_obs_norm_detail(6) + toc(t_obs_detail)
                    if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                        if( L_BENCH_GLOB ) t_obs_detail = tic()
                        weight_prev = 1. - update_frac_trail_rec
                        call vol_prev_even%mul(weight_prev)
                        call vol_prev_odd%mul(weight_prev)
                        call build%vol%mul(update_frac_trail_rec)
                        call build%vol2%mul(update_frac_trail_rec)
                        call build%vol%add(vol_prev_even)
                        call build%vol2%add(vol_prev_odd)
                        call build%vol%write(eonames(1))
                        call build%vol2%write(eonames(2))
                        if( L_BENCH_GLOB ) rt_obs_norm_detail(5) = rt_obs_norm_detail(5) + toc(t_obs_detail)
                    endif
                    if( L_BENCH_GLOB ) t_obs_detail = tic()
                    call build%vol_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                        &params%smpd_crop, wthreads=.true.)
                    call build%vol%pad_fft(build%vol_pad)
                    call build%vol_pad%expand_cmat(params%box)
                    call vol_pad2ref_pfts(build%pftc, build%vol_pad, build%eulspace, state, .true.)
                    call build%vol_pad%kill
                    call build%vol_pad%kill_expanded
                    call build%vol_odd_pad%new([params%box_croppd, params%box_croppd, params%box_croppd], &
                        &params%smpd_crop, wthreads=.true.)
                    call build%vol2%pad_fft(build%vol_odd_pad)
                    call build%vol_odd_pad%expand_cmat(params%box)
                    call vol_pad2ref_pfts(build%pftc, build%vol_odd_pad, build%eulspace, state, .false.)
                    call build%vol_odd_pad%kill
                    call build%vol_odd_pad%kill_expanded
                    if( L_BENCH_GLOB ) rt_obs_norm_detail(8) = rt_obs_norm_detail(8) + toc(t_obs_detail)
                    fsc_file = refine3D_fsc_fname(state)
                    if( allocated(fsc_state) ) deallocate(fsc_state)
                    fsc_state = file2rarr(fsc_file)
                    do iproj = 1, params%nspace
                        call obs_frcs%set_frc(iproj, fsc_state, state)
                    enddo
                    if( allocated(fsc_state) ) deallocate(fsc_state)
                    call vol_prev_even%kill
                    call vol_prev_odd%kill
                    if( L_BENCH_GLOB ) rt_normalize = rt_normalize + toc(t_normalize)
                enddo
                call obs_frcs%write(string(FRCS_FILE))
                call obs_frcs%kill
                if( L_BENCH_GLOB ) t_write = tic()
                call write_polar_refs_from_current_pftc(params, build)
                if( L_BENCH_GLOB ) rt_write = toc(t_write)
                l_refs_written = .true.
                call gridcorr_img%kill
                call vol_e%kill
                call vol_merged%kill
                call fsc_file%kill
                call fsc_txt_file%kill
                call recname%kill
                call volname%kill
                call volname_prev%kill
                call volname_prev_even%kill
                call volname_prev_odd%kill
                if( allocated(fsc) ) deallocate(fsc)
                if( allocated(fsc_state) ) deallocate(fsc_state)
                if( L_BENCH_GLOB ) rt_reduce = rt_obs_reduce_detail(1) + rt_obs_reduce_detail(2)
            case('yes')
                if( L_BENCH_GLOB ) t_reduce = tic()
                call build%pftc%polar_cavger_assemble_sums_from_parts
                if( L_BENCH_GLOB ) rt_reduce = toc(t_reduce)
                if( L_BENCH_GLOB ) t_normalize = tic()
                call build%pftc%polar_cavger_normalize_commonline_refs(build%eulspace, build%pgrpsyms, cline, update_frac_eff)
                if( L_BENCH_GLOB ) rt_normalize = toc(t_normalize)
            case default
                THROW_HARD('unsupported POLAR mode: '//trim(params%polar))
        end select
        if( L_BENCH_GLOB ) t_resolution = tic()
        call update_polar_resolution_fields(params, build)
        if( L_BENCH_GLOB ) rt_resolution = toc(t_resolution)
        if( .not. l_refs_written )then
            if( L_BENCH_GLOB ) t_write = tic()
            call build%pftc%polar_cavger_writeall(refine3D_polar_refs_fbody())
            if( L_BENCH_GLOB ) rt_write = toc(t_write)
        endif
        if( L_BENCH_GLOB ) t_cleanup = tic()
        call eorecvol_read%kill
        call build%pftc%polar_cavger_kill
        call build%pftc%kill
        call build%kill_rec_eo_tbox
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        call simple_touch('VOLASSEMBLE_FINISHED')
        if( L_BENCH_GLOB ) rt_cleanup = toc(t_cleanup)
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            rt_cmp_restore  = rt_normalize
            rt_cmp_accounted= rt_setup + rt_reduce + rt_cmp_restore + rt_resolution + rt_write + rt_cleanup
            benchfname = refine3D_volassemble_bench_fname(params%which_iter)
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a,a)')  'volassemble assembly mode           : polar=', trim(params%polar)
            write(fnr,'(a,i0)') 'volassemble nspace                  : ', params%nspace
            write(fnr,'(a,i0)') 'volassemble pftsz                   : ', params%pftsz
            write(fnr,'(a,i0)') 'volassemble nstates                 : ', params%nstates
            write(fnr,'(a,i0)') 'volassemble kfrom                   : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'volassemble kto                     : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble setup/init              : ', rt_setup
            write(fnr,'(a,t52,f9.2)') 'volassemble reduce partial inputs   : ', rt_reduce
            write(fnr,'(a,t52,f9.2)') 'volassemble restore/reconstruct refs: ', rt_cmp_restore
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking             : ', 0.
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering    : ', 0.
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata        : ', rt_resolution
            write(fnr,'(a,t52,f9.2)') 'volassemble reference handoff write : ', rt_write
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble total time              : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE DETAIL TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble obsfield volume restore : ', rt_obs_norm_detail(7)
            write(fnr,'(a,t52,f9.2)') 'volassemble obsfield polar project  : ', rt_obs_norm_detail(8)
            write(fnr,'(a,t52,f9.2)') 'volassemble previous reference read : ', rt_obs_norm_detail(4)
            write(fnr,'(a,t52,f9.2)') 'volassemble trailing blend/vols     : ', rt_obs_norm_detail(5)
            write(fnr,'(a,t52,f9.2)') 'volassemble lowres even/odd insert  : ', rt_obs_norm_detail(6)
            write(fnr,'(a,t52,f9.2)') 'volassemble FSC/FRC bookkeeping     : ', rt_obs_norm_detail(12)
            write(fnr,'(a,t52,f9.2)') 'volassemble resolution metadata     : ', rt_resolution
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', rt_cleanup
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble polar assembly setup     : ', rt_setup
            select case(trim(params%polar))
                case('obsfield')
                    write(fnr,'(a,t52,f9.2)') 'volassemble reduce obsfield rec partials: ', rt_reduce
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield read parts         : ', rt_obs_reduce_detail(1)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield sum partial recs   : ', rt_obs_reduce_detail(2)
                    write(fnr,'(a,t52,f9.2)') 'volassemble normalize obsfield refs     : ', rt_normalize
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield read prev refs     : ', rt_obs_norm_detail(4)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield trail volumes      : ', rt_obs_norm_detail(5)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield lowres ref insert  : ', rt_obs_norm_detail(6)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield volume restore     : ', rt_obs_norm_detail(7)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield polar project      : ', rt_obs_norm_detail(8)
                    write(fnr,'(a,t52,f9.2)') 'volassemble obsfield fsc/frc bookkeeping: ', rt_obs_norm_detail(12)
                case('yes')
                    write(fnr,'(a,t52,f9.2)') 'volassemble reduce polar sums           : ', rt_reduce
                    write(fnr,'(a,t52,f9.2)') 'volassemble normalize common-line refs  : ', rt_normalize
            end select
            write(fnr,'(a,t52,f9.2)') 'volassemble resolution metadata      : ', rt_resolution
            write(fnr,'(a,t52,f9.2)') 'volassemble write POLAR_REFS handoff : ', rt_write
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize         : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble polar assembly total     : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble setup/init              : ', (rt_setup/rt_tot)       * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble reduce partial inputs   : ', (rt_reduce/rt_tot)      * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble restore/reconstruct refs: ', (rt_cmp_restore/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking             : ', 0.
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering    : ', 0.
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata        : ', (rt_resolution/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble reference handoff write : ', (rt_write/rt_tot)       * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', (rt_cleanup/rt_tot)     * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble % accounted for         : ', (rt_cmp_accounted/rt_tot)* 100.
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble polar assembly setup     : ', (rt_setup/rt_tot)     * 100.
            select case(trim(params%polar))
                case('obsfield')
                    write(fnr,'(a,t52,f9.2)') 'volassemble reduce obsfield rec partials: ', (rt_reduce/rt_tot)    * 100.
                    write(fnr,'(a,t52,f9.2)') 'volassemble normalize obsfield refs     : ', (rt_normalize/rt_tot) * 100.
                case('yes')
                    write(fnr,'(a,t52,f9.2)') 'volassemble reduce polar sums           : ', (rt_reduce/rt_tot)    * 100.
                    write(fnr,'(a,t52,f9.2)') 'volassemble normalize common-line refs  : ', (rt_normalize/rt_tot) * 100.
            end select
            write(fnr,'(a,t52,f9.2)') 'volassemble resolution metadata      : ', (rt_resolution/rt_tot)* 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble write POLAR_REFS handoff : ', (rt_write/rt_tot)     * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize         : ', (rt_cleanup/rt_tot)   * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble % accounted for          : ', &
                ((rt_setup+rt_reduce+rt_normalize+rt_resolution+rt_write+rt_cleanup)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_polar_assembly

    subroutine update_polar_resolution_fields( params, build )
        type(parameters), intent(in)    :: params
        type(builder),    intent(inout) :: build
        type(string) :: fsc_file
        real,    allocatable :: fsc(:), res(:), res0143s(:)
        logical, allocatable :: has_fsc(:)
        real    :: fsc05, fsc0143
        integer :: state, iptcl, istate
        allocate(res0143s(params%nstates), source=0.)
        allocate(has_fsc(params%nstates), source=.false.)
        res = get_resarr(params%box_crop, params%smpd_crop)
        do state = 1, params%nstates
            fsc_file = refine3D_fsc_fname(state)
            if( .not. file_exists(fsc_file) ) cycle
            fsc = file2rarr(fsc_file)
            call get_resolution(fsc, res, fsc05, fsc0143)
            res0143s(state) = fsc0143
            has_fsc(state)  = .true.
            deallocate(fsc)
        end do
        if( any(has_fsc) )then
            if( params%nstates == 1 )then
                if( has_fsc(1) ) call build%spproj_field%set_all2single('res', res0143s(1))
            else
                do iptcl = 1, build%spproj_field%get_noris()
                    istate = build%spproj_field%get_state(iptcl)
                    if( istate > 0 .and. istate <= params%nstates )then
                        if( has_fsc(istate) ) call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                    endif
                end do
            endif
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        call fsc_file%kill
        if( allocated(fsc)      ) deallocate(fsc)
        if( allocated(res)      ) deallocate(res)
        if( allocated(res0143s) ) deallocate(res0143s)
        if( allocated(has_fsc)  ) deallocate(has_fsc)
    end subroutine update_polar_resolution_fields

    subroutine exec_cartesian_assembly( self, cline )
        use simple_reconstructor_eo,        only: reconstructor_eo
        use simple_gridding,                only: prep3D_inv_instrfun4mul
        use simple_matcher_smpl_and_lplims, only: set_bp_range3D
        use simple_matcher_refvol_utils,    only: read_mask_filter_reproject_refvols, write_polar_refs_from_current_pftc
        use simple_nu_filter,               only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
        &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, analyze_filtmap_neighbor_continuity
        use simple_vol_pproc_policy, only: vol_pproc_plan, plan_state_postprocess, AUTOMASK_ACTION_REGENERATE,&
                                                   &NU_MASK_SOURCE_FRESH_AUTOMASK, NU_MASK_SOURCE_EXISTING_AUTOMASK
        class(commander_cartesian_volassemble), intent(inout) :: self
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
        type(string)                  :: recname, volname, volname_prev, fsc_txt_file
        type(string)                  :: volname_prev_even, volname_prev_odd
        type(string)                  :: eonames(2), eonames_nu(2), volname_nu, benchfname, write_polar_refs_arg
        type(vol_pproc_plan) :: pp_plan
        logical, allocatable          :: l_mask(:,:,:)
        logical                       :: l_nonuniform_mode, l_write_polar_refs
        integer, allocatable          :: imat(:,:,:)
        real, allocatable             :: fsc(:), res05s(:), res0143s(:)
        real                          :: weight_prev, update_frac_trail_rec, mskrad_px
        integer                       :: part, state, iptcl, istate, find4eoavg, fnr, ldim(3), ldim_pd(3), numlen_part
        integer                       :: which_iter
        integer(timer_int_kind)       :: t_init, t_read, t_sum_reduce, t_sum_eos, t_sampl_dens_correct_eos
        integer(timer_int_kind)       :: t_sampl_dens_correct_sum, t_eoavg, t_automask, t_nonuniform
        integer(timer_int_kind)       :: t_project_refs, t_project_metadata, t_cleanup, t_tot
        integer(timer_int_kind)       :: t_prev_ref_read, t_prev_fsc_prior, t_lowres_insert, t_trailing_refs
        real(timer_int_kind)          :: rt_init, rt_read, rt_sum_reduce, rt_sum_eos, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_eoavg, rt_automask, rt_nonuniform
        real(timer_int_kind)          :: rt_project_refs, rt_project_metadata, rt_cleanup, rt_tot, rt_accounted
        real(timer_int_kind)          :: rt_prev_ref_read, rt_prev_fsc_prior, rt_lowres_insert, rt_trailing_refs
        real(timer_int_kind)          :: rt_cmp_reduce, rt_cmp_restore
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_rec_eo_tbox(params) ! reconstruction toolbox built
        call build%eorecvol%kill_exp         ! reduced memory usage
        l_write_polar_refs = .true.
        if( cline%defined('write_polar_refs') )then
            write_polar_refs_arg = cline%get_carg('write_polar_refs')
            l_write_polar_refs   = write_polar_refs_arg%to_char() /= 'no'
        endif
        numlen_part       = max(1, params%numlen)
        l_nonuniform_mode = trim(params%filt_mode).eq.'nonuniform'
        allocate(res05s(params%nstates), res0143s(params%nstates))
        res0143s = 0.
        res05s   = 0.
        call eorecvol_read%new(params, build%spproj, expand=.false.)
        if( L_BENCH_GLOB )then
            ! end of init
            rt_init = toc(t_init)
            ! initialise incremental timers before loop
            rt_read                    = 0.
            rt_sum_reduce              = 0.
            rt_sum_eos                 = 0.
            rt_sampl_dens_correct_eos  = 0.
            rt_sampl_dens_correct_sum  = 0.
            rt_eoavg                   = 0.
            rt_automask                = 0.
            rt_nonuniform              = 0.
            rt_project_refs            = 0.
            rt_project_metadata        = 0.
            rt_cleanup                 = 0.
            rt_prev_ref_read           = 0.
            rt_prev_fsc_prior          = 0.
            rt_lowres_insert           = 0.
            rt_trailing_refs           = 0.
        endif
        ! read in previous reconstruction when trail_rec==yes
        update_frac_trail_rec = 1.0
        if( params%l_trail_rec )then
            if( cline%defined('ufrac_trec') )then
                update_frac_trail_rec = params%ufrac_trec
            else
                call build%spproj%read_segment(params%oritype, params%projfile)
                update_frac_trail_rec = build%spproj%os_ptcl3D%get_update_frac()
            endif
        endif
        ! Prep for correction of the shape of the interpolator
        ldim         = build%vol%get_ldim()
        ldim_pd      = OSMPL_PAD_FAC * ldim
        gridcorr_img = prep3D_inv_instrfun4mul(ldim, ldim_pd, params%smpd_crop)
        ! assemble volumes
        do state=1,params%nstates
            call build%eorecvol%reset_all
            ! assemble volumes
            do part=1,params%nparts
                if( L_BENCH_GLOB ) t_read = tic()
                call eorecvol_read%read_eos(refine3D_partial_rec_fbody(state, part, numlen_part))
                ! sum the Fourier coefficients
                if( L_BENCH_GLOB )then
                    rt_read       = rt_read + toc(t_read)
                    t_sum_reduce  = tic()
                endif
                call build%eorecvol%sum_reduce(eorecvol_read)
                if( L_BENCH_GLOB ) rt_sum_reduce = rt_sum_reduce + toc(t_sum_reduce)
            end do
            ! correct for sampling density and estimate resolution
            recname    = refine3D_state_vol_fbody(state)
            volname    = refine3D_state_vol_fname(state)
            eonames(1) = refine3D_state_halfvol_fname(state, 'even')
            eonames(2) = refine3D_state_halfvol_fname(state, 'odd')
            if( params%l_ml_reg )then
                ! the sum is done after regularization
            else
                if( L_BENCH_GLOB ) t_sum_eos = tic()
                call build%eorecvol%sum_eos
                if( L_BENCH_GLOB ) rt_sum_eos = rt_sum_eos + toc(t_sum_eos)
            endif
            if( L_BENCH_GLOB ) t_sampl_dens_correct_eos = tic()
            if( params%l_trail_rec )then
                if( .not. cline%defined('vol'//int2str(state)) ) THROW_HARD('vol'//int2str(state)//' required in volassemble cline when trail_rec==yes')
                volname_prev      = cline%get_carg('vol'//int2str(state))
                volname_prev_even = add2fbody(volname_prev, MRC_EXT, '_even')
                volname_prev_odd  = add2fbody(volname_prev, MRC_EXT, '_odd')
                if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
                if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
                if( L_BENCH_GLOB ) t_prev_ref_read = tic()
                call vol_prev_even%read_and_crop(volname_prev_even, params%smpd, params%box_crop, params%smpd_crop)
                call vol_prev_odd%read_and_crop( volname_prev_odd,  params%smpd, params%box_crop, params%smpd_crop)
                if( L_BENCH_GLOB ) rt_prev_ref_read = rt_prev_ref_read + toc(t_prev_ref_read)
                if( allocated(fsc) ) deallocate(fsc)
                if( L_BENCH_GLOB ) t_prev_fsc_prior = tic()
                call build%eorecvol%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                if( L_BENCH_GLOB ) rt_prev_fsc_prior = rt_prev_fsc_prior + toc(t_prev_fsc_prior)
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg, fsc)
            else 
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            endif
            if( cline%defined('which_iter') )then
                fsc_txt_file = refine3D_resolution_txt_fbody(state, params%which_iter)
            else
                fsc_txt_file = refine3D_resolution_txt_fbody(state)
            endif
            call build%eorecvol%write_fsc2txt(fsc_txt_file)
            if( L_BENCH_GLOB ) rt_sampl_dens_correct_eos = rt_sampl_dens_correct_eos + toc(t_sampl_dens_correct_eos)
            if( params%l_ml_reg )then
                if( L_BENCH_GLOB ) t_sum_eos = tic()
                call build%eorecvol%sum_eos
                if( L_BENCH_GLOB ) rt_sum_eos = rt_sum_eos + toc(t_sum_eos)
            endif
            call build%eorecvol%get_res(res05s(state), res0143s(state))
            if( L_BENCH_GLOB ) t_sampl_dens_correct_sum = tic()
            call build%eorecvol%sampl_dens_correct_sum( build%vol )
            if( L_BENCH_GLOB ) rt_sampl_dens_correct_sum = rt_sampl_dens_correct_sum + toc(t_sampl_dens_correct_sum)
            ! need to put the sum back at lowres for the eo pairs
            if( L_BENCH_GLOB ) t_eoavg = tic()
            if( L_BENCH_GLOB ) t_lowres_insert = tic()
            call build%vol%fft
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
            call build%vol%copy(vol_e)
            ! merged volume in separate object to preserve even/odd in memory
            call vol_merged%copy(vol_e)
            call vol_merged%add(build%vol2)
            call vol_merged%mul(0.5)
            call vol_merged%write( volname, del_if_exists=.true. )
            call wait_for_closure( volname )
            if( L_BENCH_GLOB ) rt_lowres_insert = rt_lowres_insert + toc(t_lowres_insert)
            if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                if( L_BENCH_GLOB ) t_trailing_refs = tic()
                weight_prev = 1. - update_frac_trail_rec
                call vol_prev_even%mul(weight_prev)
                call vol_prev_odd%mul (weight_prev)
                call build%vol%mul(update_frac_trail_rec)
                call build%vol2%mul(update_frac_trail_rec)
                call build%vol%add(vol_prev_even)
                call build%vol2%add(vol_prev_odd)
                call build%vol%write(eonames(1))  ! even trailed
                call build%vol2%write(eonames(2)) ! odd trailed
                call vol_prev_even%kill
                call vol_prev_odd%kill
                if( L_BENCH_GLOB ) rt_trailing_refs = rt_trailing_refs + toc(t_trailing_refs)
            endif
            params%vols(state)      = volname
            params%vols_even(state) = eonames(1)
            params%vols_odd(state)  = eonames(2)
            if( L_BENCH_GLOB ) rt_eoavg = rt_eoavg + toc(t_eoavg)
            which_iter = 1
            if( cline%defined('which_iter') ) which_iter = params%which_iter
            call plan_state_postprocess(params, state, which_iter, l_nonuniform_mode, pp_plan)
            if( pp_plan%l_state_mask_incompatible )then
                write(logfhandle,'(A,1X,A)') '>>> Existing automask incompatible with current box/sampling, regenerating:', &
                    &pp_plan%mskfile_state%to_char()
            endif
            ! automasking after even/odd volume generation (state-aware for multi-state support)
            if( L_BENCH_GLOB ) t_automask = tic()
            if( pp_plan%automask_action == AUTOMASK_ACTION_REGENERATE )then
                call mskvol%automask3D(params, build%vol, build%vol2, pp_plan%automask_tight)
                call mskvol%write(pp_plan%mskfile_state)
            endif
            if( L_BENCH_GLOB ) rt_automask = rt_automask + toc(t_automask)
            if( l_nonuniform_mode )then
                if( L_BENCH_GLOB ) t_nonuniform = tic()
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
                if( allocated(nu_aux_even) ) deallocate(nu_aux_even)
                if( allocated(nu_aux_odd) )  deallocate(nu_aux_odd)
                if( params%l_ml_reg ) then
                    call vol_nu_base_even%new(ldim, params%smpd_crop)
                    call vol_nu_base_odd%new( ldim, params%smpd_crop)
                    call vol_nu_base_even%read(add2fbody(eonames(1), MRC_EXT, '_unfil'))
                    call vol_nu_base_odd%read( add2fbody(eonames(2), MRC_EXT, '_unfil'))
                    allocate(nu_aux_even(1), nu_aux_odd(1))
                    call nu_aux_even(1)%copy(build%vol)
                    call nu_aux_odd(1)%copy(build%vol2)
                    call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, l_mask, [res0143s(state)], &
                        &nu_aux_even, nu_aux_odd)
                else
                    ! build%vol/build%vol2 hold the current even/odd pair in memory.
                    call setup_nu_dmats(build%vol, build%vol2, l_mask, [real ::])
                end if
                call optimize_nu_cutoff_finds()
                call nu_filter_vols(vol_even_nu, vol_odd_nu)
                if( allocated(nu_aux_even) ) then
                    call print_nu_filtmap_lowpass_stats(l_mask, aux_resolutions=[res0143s(state)])
                else
                    call print_nu_filtmap_lowpass_stats(l_mask)
                endif
                call analyze_filtmap_neighbor_continuity(l_mask)
                eonames_nu(1) = add2fbody(eonames(1), MRC_EXT, NUFILT_SUFFIX)
                eonames_nu(2) = add2fbody(eonames(2), MRC_EXT, NUFILT_SUFFIX)
                volname_nu    = add2fbody(volname,    MRC_EXT, NUFILT_SUFFIX)
                call vol_even_nu%write(eonames_nu(1), del_if_exists=.true.)
                call vol_odd_nu%write(eonames_nu(2), del_if_exists=.true.)
                call vol_even_nu%add(vol_odd_nu)
                call vol_even_nu%mul(0.5)
                call vol_even_nu%write(volname_nu, del_if_exists=.true.)
                call wait_for_closure(volname_nu)
                call vol_even_nu%kill
                call vol_odd_nu%kill
                call vol_nu_base_even%kill
                call vol_nu_base_odd%kill
                if( allocated(nu_aux_even) ) then
                    call nu_aux_even(1)%kill
                    deallocate(nu_aux_even)
                end if
                if( allocated(nu_aux_odd) ) then
                    call nu_aux_odd(1)%kill
                    deallocate(nu_aux_odd)
                end if
                call vol_msk%kill
                deallocate(l_mask)
                call cleanup_nu_filter()
                if( L_BENCH_GLOB ) rt_nonuniform = rt_nonuniform + toc(t_nonuniform)
            endif
            call recname%kill
            call volname%kill
        end do
        ! Update per-particle FSC(0.143) resolution using the values computed in assembly.
        if( L_BENCH_GLOB ) t_project_metadata = tic()
        if( params%nstates == 1 )then
            call build%spproj_field%set_all2single('res', res0143s(1))
        else
            do iptcl = 1, build%spproj_field%get_noris()
                istate = build%spproj_field%get_state(iptcl)
                if( istate > 0 .and. istate <= params%nstates )then
                    call build%spproj_field%set(iptcl, 'res', res0143s(istate))
                endif
            end do
        endif
        if( L_BENCH_GLOB ) rt_project_metadata = rt_project_metadata + toc(t_project_metadata)
        ! Cartesian refinement still matches in polar central-section space.
        ! Refinement assembly refreshes the projected reference sections for the
        ! next matcher/probability-table pass. Standalone/final reconstruct3D
        ! calls assemble only volumes and do not own matcher handoff state.
        if( l_write_polar_refs )then
            if( L_BENCH_GLOB ) t_project_refs = tic()
            call set_bp_range3D(params, build, cline)
            call read_mask_filter_reproject_refvols(params, build, cline, 1)
            call write_polar_refs_from_current_pftc(params, build)
            call build%pftc%polar_cavger_kill
            call build%pftc%kill
            if( L_BENCH_GLOB ) rt_project_refs = rt_project_refs + toc(t_project_refs)
        endif
        if( L_BENCH_GLOB ) t_project_metadata = tic()
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        if( L_BENCH_GLOB ) rt_project_metadata = rt_project_metadata + toc(t_project_metadata)
        ! destruct
        if( L_BENCH_GLOB ) t_cleanup = tic()
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
        if( allocated(nu_aux_even) ) then
            call nu_aux_even(1)%kill
            deallocate(nu_aux_even)
        end if
        if( allocated(nu_aux_odd) ) then
            call nu_aux_odd(1)%kill
            deallocate(nu_aux_odd)
        end if
        if( allocated(l_mask) ) deallocate(l_mask)
        if( allocated(imat) ) deallocate(imat)
        call cleanup_nu_filter()
        call state_mask_bin%kill_bimg
        call mskvol%kill_bimg
        call vol_e%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        call simple_touch('VOLASSEMBLE_FINISHED')
        if( L_BENCH_GLOB ) rt_cleanup = toc(t_cleanup)
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            rt_cmp_reduce  = rt_read + rt_sum_reduce
            rt_cmp_restore = rt_sum_eos + rt_sampl_dens_correct_eos + rt_sampl_dens_correct_sum + rt_eoavg
            benchfname = refine3D_volassemble_bench_fname(params%which_iter)
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** BENCHMARK CONTEXT ***'
            write(fnr,'(a)')    'volassemble assembly mode           : polar=no'
            write(fnr,'(a,i0)') 'volassemble nspace                  : ', params%nspace
            write(fnr,'(a,i0)') 'volassemble nstates                 : ', params%nstates
            write(fnr,'(a,i0)') 'volassemble kfrom                   : ', params%kfromto(1)
            write(fnr,'(a,i0)') 'volassemble kto                     : ', params%kfromto(2)
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble setup/init              : ', rt_init
            write(fnr,'(a,t52,f9.2)') 'volassemble reduce partial inputs   : ', rt_cmp_reduce
            write(fnr,'(a,t52,f9.2)') 'volassemble restore/reconstruct refs: ', rt_cmp_restore
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking             : ', rt_automask
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering    : ', rt_nonuniform
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata        : ', rt_project_metadata
            write(fnr,'(a,t52,f9.2)') 'volassemble reference handoff write : ', rt_project_refs
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble total time              : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE DETAIL TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble previous FSC/prior      : ', rt_prev_fsc_prior
            write(fnr,'(a,t52,f9.2)') 'volassemble previous reference read : ', rt_prev_ref_read
            write(fnr,'(a,t52,f9.2)') 'volassemble trailing blend/refs     : ', rt_trailing_refs
            write(fnr,'(a,t52,f9.2)') 'volassemble lowres even/odd insert  : ', rt_lowres_insert
            write(fnr,'(a,t52,f9.2)') 'volassemble FSC/FRC bookkeeping     : ', rt_sampl_dens_correct_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble resolution metadata     : ', rt_project_metadata
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', rt_cleanup
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble initialisation           : ', rt_init
            write(fnr,'(a,t52,f9.2)') 'volassemble reading of volumes (I/O) : ', rt_read
            write(fnr,'(a,t52,f9.2)') 'volassemble summing partial volumes  : ', rt_sum_reduce
            write(fnr,'(a,t52,f9.2)') 'volassemble sum of eo-paris          : ', rt_sum_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble gridding correction (eos): ', rt_sampl_dens_correct_eos
            write(fnr,'(a,t52,f9.2)') 'volassemble gridding correction (sum): ', rt_sampl_dens_correct_sum
            write(fnr,'(a,t52,f9.2)') 'volassemble averaging eo-pairs       : ', rt_eoavg
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking              : ', rt_automask
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering     : ', rt_nonuniform
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata         : ', rt_project_metadata
            write(fnr,'(a,t52,f9.2)') 'volassemble polar ref projection     : ', rt_project_refs
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize         : ', rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble total time               : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** COMPARABLE RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble setup/init              : ', (rt_init/rt_tot)         * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble reduce partial inputs   : ', (rt_cmp_reduce/rt_tot)   * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble restore/reconstruct refs: ', (rt_cmp_restore/rt_tot)  * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking             : ', (rt_automask/rt_tot)     * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering    : ', (rt_nonuniform/rt_tot)   * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata        : ', (rt_project_metadata/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble reference handoff write : ', (rt_project_refs/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize        : ', (rt_cleanup/rt_tot)     * 100.
            rt_accounted = rt_init + rt_cmp_reduce + rt_cmp_restore + rt_automask + rt_nonuniform + &
                &rt_project_metadata + rt_project_refs + rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble % accounted for         : ', (rt_accounted/rt_tot)    * 100.
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,t52,f9.2)') 'volassemble initialisation           : ', (rt_init/rt_tot)                   * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble reading of volumes (I/O) : ', (rt_read/rt_tot)                   * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble summing partial volumes  : ', (rt_sum_reduce/rt_tot)             * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble sum of eo-paris          : ', (rt_sum_eos/rt_tot)                * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble gridding correction (eos): ', (rt_sampl_dens_correct_eos/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble gridding correction (sum): ', (rt_sampl_dens_correct_sum/rt_tot) * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble averaging eo-pairs       : ', (rt_eoavg/rt_tot)                  * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble automasking              : ', (rt_automask/rt_tot)               * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble nonuniform filtering     : ', (rt_nonuniform/rt_tot)             * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble project metadata         : ', (rt_project_metadata/rt_tot)       * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble polar ref projection     : ', (rt_project_refs/rt_tot)          * 100.
            write(fnr,'(a,t52,f9.2)') 'volassemble cleanup/finalize         : ', (rt_cleanup/rt_tot)               * 100.
            rt_accounted = rt_init + rt_read + rt_sum_reduce + rt_sum_eos + rt_sampl_dens_correct_eos + &
                &rt_sampl_dens_correct_sum + rt_eoavg + rt_automask + rt_nonuniform + rt_project_metadata + &
                &rt_project_refs + rt_cleanup
            write(fnr,'(a,t52,f9.2)') 'volassemble % accounted for          : ', (rt_accounted/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_cartesian_assembly

end module simple_commanders_rec_distr
