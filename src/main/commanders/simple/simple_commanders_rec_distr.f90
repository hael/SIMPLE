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
        class(commander_volassemble), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(reconstructor_eo)        :: eorecvol_read
        type(image)                   :: vol_prev_even, vol_prev_odd, gridcorr_img
        type(string)                  :: recname, volname, volname_prev, fsc_txt_file
        type(string)                  :: volname_prev_even, volname_prev_odd, str_state, str_iter
        type(string)                  :: eonames(2), benchfname
        real, allocatable             :: fsc(:), res05s(:), res0143s(:)
        real                          :: weight_prev, update_frac_trail_rec
        integer                       :: part, state, find4eoavg, fnr, ldim(3), ldim_pd(3)
        integer(timer_int_kind)       :: t_init, t_read, t_sum_reduce, t_sum_eos, t_sampl_dens_correct_eos
        integer(timer_int_kind)       :: t_sampl_dens_correct_sum, t_eoavg, t_tot
        real(timer_int_kind)          :: rt_init, rt_read, rt_sum_reduce, rt_sum_eos, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_eoavg, rt_tot
        if( L_BENCH_GLOB )then
            t_init = tic()
            t_tot  = t_init
        endif
        call build%init_params_and_build_general_tbox(cline,params)
        call build%build_rec_eo_tbox(params) ! reconstruction toolbox built
        call build%eorecvol%kill_exp         ! reduced memory usage
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
                call eorecvol_read%read_eos(string(VOL_FBODY)//int2str_pad(state,2)//'_part'//int2str_pad(part,params%numlen))
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
            call build%vol2%zero_and_unflag_ft
            call build%vol2%read(eonames(2))
            call build%vol2%fft()
            call build%vol2%insert_lowres(build%vol, find4eoavg)
            call build%vol2%ifft()
            call build%vol2%mul(gridcorr_img)
            call build%vol2%write(eonames(2), del_if_exists=.true.)
            ! merged volume
            call build%vol%ifft
            call build%vol%mul(gridcorr_img)
            call build%vol%write( volname, del_if_exists=.true. )
            call wait_for_closure( volname )
            if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                call build%vol%read(eonames(1))  ! even current
                call build%vol2%read(eonames(2)) ! odd current
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
            write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
            &((rt_init+rt_read+rt_sum_reduce+rt_sum_eos+rt_sampl_dens_correct_eos+rt_sampl_dens_correct_sum+rt_eoavg)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_volassemble

end module simple_commanders_rec_distr
