! concrete commander: 3D reconstruction routines
module simple_commander_rec
include 'simple_lib.f08'
use simple_parameters,          only: parameters
use simple_builder,             only: builder
use simple_cmdline,             only: cmdline
use simple_commander_base,      only: commander_base
use simple_projection_frcs,     only: projection_frcs
use simple_strategy2D3D_common, only: gen_projection_frcs
implicit none

public :: reconstruct3D_commander
public :: volassemble_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: reconstruct3D_commander
  contains
    procedure :: execute      => exec_reconstruct3D
end type reconstruct3D_commander
type, extends(commander_base) :: volassemble_commander
  contains
    procedure :: execute      => exec_volassemble
end type volassemble_commander

contains

    !> for reconstructing volumes from image stacks and their estimated orientations
    subroutine exec_reconstruct3D( self, cline )
        use simple_rec_master, only: exec_rec, exec_rec_soft
        class(reconstruct3D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_rec_eo_tbox(params)
        if( params%l_rec_soft )then
            call build%build_strategy3D_tbox(params)
            call exec_rec_soft(cline, 1) ! which_iter = 1
        else
            call exec_rec
        endif
        ! end gracefully
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D

    !> for assembling even/odd volumes generated with distributed execution
    subroutine exec_volassemble( self, cline )
        use simple_reconstructor_eo, only: reconstructor_eo
        use simple_filterer,         only: gen_anisotropic_optlp
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(reconstructor_eo)        :: eorecvol_read
        character(len=:), allocatable :: finished_fname, recname, volname
        character(len=LONGSTRLEN)     :: eonames(2), resmskname, benchfname
        real, allocatable             :: res05s(:), res0143s(:)
        real                          :: res
        integer                       :: part, s, n, ss, state, find4eoavg, fnr
        logical, parameter            :: L_BENCH = .false.
        integer(timer_int_kind)       :: t_init, t_assemble, t_sum_eos, t_sampl_dens_correct_eos
        integer(timer_int_kind)       :: t_gen_projection_frcs, t_gen_anisotropic_optlp
        integer(timer_int_kind)       :: t_sampl_dens_correct_sum, t_eoavg, t_tot
        real(timer_int_kind)          :: rt_init, rt_assemble, rt_sum_eos, rt_sampl_dens_correct_eos
        real(timer_int_kind)          :: rt_gen_projection_frcs, rt_gen_anisotropic_optlp
        real(timer_int_kind)          :: rt_sampl_dens_correct_sum, rt_eoavg, rt_tot
        if( L_BENCH )then
            t_init = tic()
            t_tot  = t_init
        endif
        call build%init_params_and_build_general_tbox(cline,params)
        call build%build_rec_eo_tbox(params) ! reconstruction toolbox built
        call build%eorecvol%kill_exp         ! reduced meory usage
        allocate(res05s(params%nstates), res0143s(params%nstates), stat=alloc_stat)
        if(alloc_stat.ne.0)call allocchk("In: simple_eo_volassemble res05s res0143s",alloc_stat)
        res0143s = 0.
        res05s   = 0.
        ! rebuild build%vol according to box size (because it is otherwise boxmatch)
        call build%vol%new([params%box,params%box,params%box], params%smpd)
        call eorecvol_read%new( build%spproj)
        call eorecvol_read%kill_exp ! reduced memory usage
        n = params%nstates*params%nparts
        if( L_BENCH )then
            ! end of init
            rt_init = toc(t_init)
            ! initialise incremental timers before loop
            rt_assemble                = 0.
            rt_sum_eos                 = 0.
            rt_sampl_dens_correct_eos  = 0.
            rt_gen_projection_frcs     = 0.
            rt_gen_anisotropic_optlp   = 0.
            rt_sampl_dens_correct_sum  = 0.
            rt_eoavg                   = 0.
        endif
        do ss=1,params%nstates
            if( cline%defined('state') )then
                s     = 1             ! index in reconstruct3D
                state = params%state  ! actual state
            else
                s     = ss
                state = ss
            endif
            if( L_BENCH ) t_assemble = tic()
            call build%eorecvol%reset_all
            ! assemble volumes
            do part=1,params%nparts
                call eorecvol_read%read_eos(trim(VOL_FBODY)//int2str_pad(state,2)//'_part'//int2str_pad(part,params%numlen))
                ! sum the Fourier coefficients
                call build%eorecvol%sum_reduce(eorecvol_read)
            end do
            if( L_BENCH ) rt_assemble = rt_assemble + toc(t_assemble)
            ! correct for sampling density and estimate resolution
            allocate(recname, source=trim(VOL_FBODY)//int2str_pad(state,2))
            allocate(volname, source=recname//params%ext)
            eonames(1) = trim(recname)//'_even'//params%ext
            eonames(2) = trim(recname)//'_odd'//params%ext
            resmskname = params%mskfile
            if( L_BENCH ) t_sum_eos = tic()
            call build%eorecvol%sum_eos
            if( L_BENCH )then
                rt_sum_eos               = rt_sum_eos + toc(t_sum_eos)
                t_sampl_dens_correct_eos = tic()
            endif
            call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            if( L_BENCH )then
                rt_sampl_dens_correct_eos = rt_sampl_dens_correct_eos + toc(t_sampl_dens_correct_eos)
                t_gen_projection_frcs     = tic()
            endif
            call gen_projection_frcs( cline, eonames(1), eonames(2), resmskname, s, build%projfrcs)
            if( L_BENCH ) rt_gen_projection_frcs = rt_gen_projection_frcs + toc(t_gen_projection_frcs)
            call build%projfrcs%write(FRCS_FILE)
            if( L_BENCH ) t_gen_anisotropic_optlp = tic()
            call gen_anisotropic_optlp(build%vol2, build%projfrcs, build%eulspace_red, s, &
                &params%pgrp, params%hpind_fsc, params%l_phaseplate)
            if( L_BENCH ) rt_gen_anisotropic_optlp = rt_gen_anisotropic_optlp + toc(t_gen_anisotropic_optlp)
            call build%vol2%write('aniso_optlp_state'//int2str_pad(state,2)//params%ext)
            call build%eorecvol%get_res(res05s(s), res0143s(s))
            if( L_BENCH ) t_sampl_dens_correct_sum = tic()
            call build%eorecvol%sampl_dens_correct_sum( build%vol )
            if( L_BENCH ) rt_sampl_dens_correct_sum = rt_sampl_dens_correct_sum + toc(t_sampl_dens_correct_sum)
            call build%vol%write( volname, del_if_exists=.true. )
            call wait_for_closure( volname )
            ! need to put the sum back at lowres for the eo pairs
            if( L_BENCH ) t_eoavg = tic()
            call build%vol%fft()
            call build%vol2%zero_and_unflag_ft
            call build%vol2%read(eonames(1))
            call build%vol2%fft()
            call build%vol2%insert_lowres(build%vol, find4eoavg)
            call build%vol2%ifft()
            call build%vol2%write(eonames(1), del_if_exists=.true.)
            call build%vol2%zero_and_unflag_ft
            call build%vol2%read(eonames(2))
            call build%vol2%fft()
            call build%vol2%insert_lowres(build%vol, find4eoavg)
            call build%vol2%ifft()
            call build%vol2%write(eonames(2), del_if_exists=.true.)
            if( L_BENCH ) rt_eoavg = rt_eoavg + toc(t_eoavg)
            deallocate(recname, volname)
            if( cline%defined('state') )exit
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        params%lp = max( params%lpstop,res )
        write(logfhandle,'(a,1x,F6.2)') '>>> LOW-PASS LIMIT:', params%lp
        call eorecvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        if( cline%defined('state') )then
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2))
        else
            allocate( finished_fname, source='VOLASSEMBLE_FINISHED' )
        endif
        call simple_touch( finished_fname , errmsg='In: commander_rec::volassemble')
        if( L_BENCH )then
            rt_tot     = toc(t_tot)
            benchfname = 'volassemble_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
            write(fnr,'(a)') '*** TIMINGS (s) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation           : ', rt_init
            write(fnr,'(a,1x,f9.2)') 'assemble of volumes (I/O): ', rt_assemble
            write(fnr,'(a,1x,f9.2)') 'sum of eo-paris          : ', rt_sum_eos
            write(fnr,'(a,1x,f9.2)') 'gridding correction (eos): ', rt_sampl_dens_correct_eos
            write(fnr,'(a,1x,f9.2)') 'projection FRCs          : ', rt_gen_projection_frcs
            write(fnr,'(a,1x,f9.2)') 'anisotropic filter       : ', rt_gen_anisotropic_optlp
            write(fnr,'(a,1x,f9.2)') 'gridding correction (sum): ', rt_sampl_dens_correct_sum
            write(fnr,'(a,1x,f9.2)') 'averaging eo-pairs       : ', rt_eoavg
            write(fnr,'(a,1x,f9.2)') 'total time               : ', rt_tot
            write(fnr,'(a)') ''
            write(fnr,'(a)') '*** RELATIVE TIMINGS (%) ***'
            write(fnr,'(a,1x,f9.2)') 'initialisation           : ', (rt_init/rt_tot)                   * 100.
            write(fnr,'(a,1x,f9.2)') 'assemble of volumes (I/O): ', (rt_assemble/rt_tot)               * 100.
            write(fnr,'(a,1x,f9.2)') 'sum of eo-paris          : ', (rt_sum_eos/rt_tot)                * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding correction (eos): ', (rt_sampl_dens_correct_eos/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'projection FRCs          : ', (rt_gen_projection_frcs/rt_tot)    * 100.
            write(fnr,'(a,1x,f9.2)') 'anisotropic filter       : ', (rt_gen_anisotropic_optlp/rt_tot)  * 100.
            write(fnr,'(a,1x,f9.2)') 'gridding correction (sum): ', (rt_sampl_dens_correct_sum/rt_tot) * 100.
            write(fnr,'(a,1x,f9.2)') 'averaging eo-pairs       : ', (rt_eoavg/rt_tot)                  * 100.
            write(fnr,'(a,1x,f9.2)') '% accounted for          : ',&
            &((rt_init+rt_assemble+rt_sum_eos+rt_sampl_dens_correct_eos+rt_gen_projection_frcs+&
            &rt_gen_anisotropic_optlp+rt_sampl_dens_correct_sum+rt_eoavg)/rt_tot) * 100.
            call fclose(fnr)
        endif
    end subroutine exec_volassemble

end module simple_commander_rec
