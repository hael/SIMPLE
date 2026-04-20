module simple_commanders_rec_distr
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_volassemble
  contains
    procedure :: execute      => exec_volassemble
end type commander_volassemble

contains

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor_eo, only: reconstructor_eo
        use simple_gridding,         only: prep3D_inv_instrfun4mul
        use simple_nu_filter,        only: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, &
                                         &cleanup_nu_filter, print_nu_filtmap_lowpass_stats
        use simple_volume_postprocess_policy, only: volume_postprocess_plan, plan_state_postprocess
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
        type(string)                  :: recname, volname, volname_prev, fsc_txt_file
        type(string)                  :: volname_prev_even, volname_prev_odd, str_state, str_iter
        type(string)                  :: eonames(2), eonames_nu(2), volname_nu, benchfname
        type(volume_postprocess_plan) :: pp_plan
        logical, allocatable          :: l_mask(:,:,:)
        logical                       :: l_nonuniform_mode
        integer, allocatable          :: imat(:,:,:)
        real, allocatable             :: fsc(:), res05s(:), res0143s(:)
        real                          :: weight_prev, update_frac_trail_rec, mskrad_px
        integer                       :: part, state, find4eoavg, fnr, ldim(3), ldim_pd(3), numlen_part, which_iter
        integer(timer_int_kind)       :: t_init, t_read, t_sum_reduce, t_sum_eos, t_sampl_dens_correct_eos
        integer(timer_int_kind)       :: t_sampl_dens_correct_sum, t_eoavg, t_automask, t_nonuniform, t_tot
        real(timer_int_kind)          :: rt_init, rt_read, rt_sum_reduce, rt_sum_eos, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_eoavg, rt_automask, rt_nonuniform, rt_tot
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call build%init_params_and_build_general_tbox(cline,params)
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
            if( pp_plan%regenerate_automask )then
                call mskvol%automask3D(params, build%vol, build%vol2, trim(params%automsk).eq.'tight')
                call mskvol%write(pp_plan%mskfile_state)
            endif
            if( L_BENCH_GLOB ) rt_automask = rt_automask + toc(t_automask)
            if( l_nonuniform_mode )then
                if( L_BENCH_GLOB ) t_nonuniform = tic()
                if( allocated(l_mask) ) deallocate(l_mask)
                if( pp_plan%use_state_mask_for_nonuniform )then
                    if( pp_plan%regenerate_automask )then
                        call mskvol%set_imat
                        call mskvol%get_imat(imat)
                    else
                        call state_mask_bin%new_bimg(ldim, params%smpd_crop)
                        call state_mask_bin%read_bimg(pp_plan%mskfile_state)
                        call state_mask_bin%get_imat(imat)
                        call state_mask_bin%kill_bimg
                    endif
                    if( allocated(imat) )then
                        allocate(l_mask(ldim(1),ldim(2),ldim(3)))
                        l_mask = imat > 0
                        deallocate(imat)
                    endif
                endif
                if( .not. allocated(l_mask) )then
                    mskrad_px = 0.5 * params%mskdiam / params%smpd_crop
                    call vol_msk%disc(ldim, params%smpd_crop, mskrad_px, l_mask)
                endif
                if( allocated(nu_aux_even) ) deallocate(nu_aux_even)
                if( allocated(nu_aux_odd) )  deallocate(nu_aux_odd)
                if( params%l_ml_reg ) then
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
                call print_nu_filtmap_lowpass_stats(l_mask)
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
            write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
            &((rt_init+rt_read+rt_sum_reduce+rt_sum_eos+rt_sampl_dens_correct_eos+rt_sampl_dens_correct_sum+rt_eoavg+rt_automask+rt_nonuniform)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_volassemble

end module simple_commanders_rec_distr
