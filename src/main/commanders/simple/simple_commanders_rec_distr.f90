module simple_commanders_rec_distr
use simple_commanders_api
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

contains

    subroutine exec_polar_assembly( self, cline )
        use simple_matcher_smpl_and_lplims, only: set_bp_range3D
        class(commander_polar_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(string) :: benchfname
        integer(timer_int_kind) :: t_tot
        real(timer_int_kind)    :: rt_tot
        real                   :: update_frac_eff
        integer :: fnr, nrefs
        call build%init_params_and_build_general_tbox(cline, params)
        if( L_BENCH_GLOB ) t_tot = tic()
        ! Matchers write partition-local Cartesian partial reconstructions
        ! (polar=no) or polar partial sums (polar=yes|obsfield), while
        ! the assembly commander owns shared-memory reduction and reference update.
        call set_bp_range3D(params, build, cline)
        nrefs = params%nspace * params%nstates
        call build%pftc%new(params, nrefs, [1,1], params%kfromto)
        params%refs = string(CAVGS_ITER_FBODY)//int2str_pad(params%which_iter,3)//params%ext%to_char()
        call build%pftc%polar_cavger_new(.true., nrefs=nrefs)
        call build%pftc%polar_cavger_calc_pops(build%spproj)
        call build%pftc%polar_cavger_assemble_sums_from_parts
        update_frac_eff = params%update_frac
        if( params%l_trail_rec .and. (.not. cline%defined('ufrac_trec')) )then
            if( build%spproj_field%has_been_sampled() )then
                update_frac_eff = build%spproj_field%get_update_frac()
            endif
        endif
        select case(trim(params%polar))
            case('obsfield')
                call build%pftc%polar_cavger_merge_eos_and_norm_obsfield(build%eulspace, cline, update_frac_eff)
            case('yes')
                call build%pftc%polar_cavger_merge_eos_and_norm(build%eulspace, build%pgrpsyms, cline, update_frac_eff)
            case default
                THROW_HARD('unsupported POLAR mode: '//trim(params%polar))
        end select
        call update_polar_resolution_fields(params, build)
        call build%pftc%polar_cavger_writeall(string(POLAR_REFS_FBODY))
        call build%pftc%polar_cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'VOLASSEMBLE_BENCH_ITER'//int2str_pad(params%which_iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'polar reference assembly : ', rt_tot
            call fclose(fnr)
        endif
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        call simple_touch('VOLASSEMBLE_FINISHED')
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
            fsc_file = string(FSC_FBODY)//int2str_pad(state,2)//BIN_EXT
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
        use simple_reconstructor_eo, only: reconstructor_eo
        use simple_gridding,         only: prep3D_inv_instrfun4mul
        use simple_matcher_smpl_and_lplims, only: set_bp_range3D
        use simple_matcher_refvol_utils, only: read_mask_filter_reproject_refvols, &
            &write_polar_refs_from_current_pftc
        use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
                                         &cleanup_nu_filter, print_nu_filtmap_lowpass_stats, &
                                         &analyze_filtmap_neighbor_continuity
        use simple_volume_postprocess_policy, only: volume_postprocess_plan, plan_state_postprocess, &
                                                   &AUTOMASK_ACTION_REGENERATE, &
                                                   &NU_MASK_SOURCE_FRESH_AUTOMASK, &
                                                   &NU_MASK_SOURCE_EXISTING_AUTOMASK
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
        type(string)                  :: volname_prev_even, volname_prev_odd, str_state, str_iter
        type(string)                  :: eonames(2), eonames_nu(2), volname_nu, benchfname
        type(volume_postprocess_plan) :: pp_plan
        logical, allocatable          :: l_mask(:,:,:)
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: imat(:,:,:)
        real, allocatable             :: fsc(:), res05s(:), res0143s(:)
        real                          :: weight_prev, update_frac_trail_rec, mskrad_px
        integer                       :: part, state, iptcl, istate, find4eoavg, fnr, ldim(3), ldim_pd(3), numlen_part
        integer                       :: which_iter
        integer(timer_int_kind)       :: t_init, t_read, t_sum_reduce, t_sum_eos, t_sampl_dens_correct_eos
        integer(timer_int_kind)       :: t_sampl_dens_correct_sum, t_eoavg, t_automask, t_nonuniform
        integer(timer_int_kind)       :: t_project_refs, t_tot
        real(timer_int_kind)          :: rt_init, rt_read, rt_sum_reduce, rt_sum_eos, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_eoavg, rt_automask, rt_nonuniform
        real(timer_int_kind)          :: rt_project_refs, rt_tot, rt_accounted
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_rec_eo_tbox(params) ! reconstruction toolbox built
        call build%eorecvol%kill_exp         ! reduced memory usage
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
                call eorecvol_read%read_eos(string(VOL_FBODY)//int2str_pad(state,2)//'_part'//int2str_pad(part,numlen_part))
                ! sum the Fourier coefficients
                if( L_BENCH_GLOB )then
                    rt_read       = rt_read + toc(t_read)
                    t_sum_reduce  = tic()
                endif
                call build%eorecvol%sum_reduce(eorecvol_read)
                if( L_BENCH_GLOB ) rt_sum_reduce = rt_sum_reduce + toc(t_sum_reduce)
            end do
            ! correct for sampling density and estimate resolution
            recname    = VOL_FBODY//int2str_pad(state,2)
            volname    = recname//params%ext
            eonames(1) = recname//'_even'//params%ext%to_char()
            eonames(2) = recname//'_odd'//params%ext%to_char()
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
                volname_prev_even = add2fbody(volname_prev, params%ext, '_even')
                volname_prev_odd  = add2fbody(volname_prev, params%ext, '_odd')
                if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//volname_prev_even%to_char()//' does not exist!')
                if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//volname_prev_odd%to_char()//' does not exist!')
                call vol_prev_even%read_and_crop(volname_prev_even, params%smpd, params%box_crop, params%smpd_crop)
                call vol_prev_odd%read_and_crop( volname_prev_odd,  params%smpd, params%box_crop, params%smpd_crop)
                if( allocated(fsc) ) deallocate(fsc)
                call build%eorecvol%calc_fsc4sampl_dens_correct(vol_prev_even, vol_prev_odd, fsc)
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg, fsc)
            else 
                call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            endif
            str_state = int2str_pad(state,2)
            if( cline%defined('which_iter') )then
                str_iter     = int2str_pad(params%which_iter,3)
                fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()//'_ITER'//str_iter%to_char()
            else
                fsc_txt_file = 'RESOLUTION_STATE'//str_state%to_char()
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
            if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
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
                    call vol_nu_base_even%read(add2fbody(eonames(1), params%ext, '_unfil'))
                    call vol_nu_base_odd%read( add2fbody(eonames(2), params%ext, '_unfil'))
                    allocate(nu_aux_even(1), nu_aux_odd(1))
                    call nu_aux_even(1)%copy(build%vol)
                    call nu_aux_odd(1)%copy(build%vol2)
                    call setup_nu_dmats(vol_nu_base_even, vol_nu_base_odd, l_mask, nu_aux_even, nu_aux_odd)
                else
                    ! build%vol/build%vol2 hold the current even/odd pair in memory.
                    call setup_nu_dmats(build%vol, build%vol2, l_mask)
                end if
                call optimize_nu_cutoff_finds()
                call nu_filter_vols(vol_even_nu, vol_odd_nu)
                if( allocated(nu_aux_even) ) then
                    call print_nu_filtmap_lowpass_stats(l_mask, aux_resolutions=[res0143s(state)])
                else
                    call print_nu_filtmap_lowpass_stats(l_mask)
                endif
                call analyze_filtmap_neighbor_continuity(l_mask)
                eonames_nu(1) = add2fbody(eonames(1), params%ext, NUFILT_SUFFIX)
                eonames_nu(2) = add2fbody(eonames(2), params%ext, NUFILT_SUFFIX)
                volname_nu    = add2fbody(volname,    params%ext, NUFILT_SUFFIX)
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
        ! Cartesian refinement still matches in polar central-section space.
        ! Therefore Cartesian assembly always refreshes the projected reference
        ! sections consumed by the next matcher/probability-table pass.
        if( L_BENCH_GLOB ) t_project_refs = tic()
        call set_bp_range3D(params, build, cline)
        call read_mask_filter_reproject_refvols(params, build, cline, 1)
        call write_polar_refs_from_current_pftc(params, build)
        call build%pftc%polar_cavger_kill
        call build%pftc%kill
        if( L_BENCH_GLOB ) rt_project_refs = rt_project_refs + toc(t_project_refs)
        call build%spproj%write_segment_inside(params%oritype, params%projfile)
        ! destruct
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
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'VOLASSEMBLE_BENCH_ITER'//int2str_pad(params%which_iter,3)//'.txt'
            call fopen(fnr, FILE=benchfname, STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation           : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'reading of volumes (I/O) : ', rt_read
            write(fnr,'(a,1x,f9.2)') 'summing partial volumes  : ', rt_sum_reduce
            write(fnr,'(a,1x,f9.2)') 'sum of eo-paris          : ', rt_sum_eos
            write(fnr,'(a,1x,f9.2)') 'gridding correction (eos): ', rt_sampl_dens_correct_eos
            write(fnr,'(a,1x,f9.2)') 'gridding correction (sum): ', rt_sampl_dens_correct_sum
            write(fnr,'(a,1x,f9.2)') 'averaging eo-pairs       : ', rt_eoavg
            write(fnr,'(a,1x,f9.2)') 'automasking              : ', rt_automask
            write(fnr,'(a,1x,f9.2)') 'nonuniform filtering     : ', rt_nonuniform
            write(fnr,'(a,1x,f9.2)') 'polar ref projection     : ', rt_project_refs
            write(fnr,'(a,1x,f9.2)') 'total time               : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation           : ', (rt_init/rt_tot)                   * 100.
            write(fnr,'(a,1x,f9.2)') 'reading of volumes (I/O) : ', (rt_read/rt_tot)                   * 100.
            write(fnr,'(a,1x,f9.2)') 'summing partial volumes  : ', (rt_sum_reduce/rt_tot)             * 100.
            write(fnr,'(a,1x,f9.2)') 'sum of eo-paris          : ', (rt_sum_eos/rt_tot)                * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding correction (eos): ', (rt_sampl_dens_correct_eos/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding correction (sum): ', (rt_sampl_dens_correct_sum/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'averaging eo-pairs       : ', (rt_eoavg/rt_tot)                  * 100.
            write(fnr,'(a,1x,f9.2)') 'automasking              : ', (rt_automask/rt_tot)               * 100.
            write(fnr,'(a,1x,f9.2)') 'nonuniform filtering     : ', (rt_nonuniform/rt_tot)             * 100.
            write(fnr,'(a,1x,f9.2)') 'polar ref projection     : ', (rt_project_refs/rt_tot)          * 100.
            rt_accounted = rt_init + rt_read + rt_sum_reduce + rt_sum_eos + rt_sampl_dens_correct_eos + &
                &rt_sampl_dens_correct_sum + rt_eoavg + rt_automask + rt_nonuniform + rt_project_refs
            write(fnr,'(a,1x,f9.2)') '% accounted for          : ', (rt_accounted/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_cartesian_assembly

end module simple_commanders_rec_distr
