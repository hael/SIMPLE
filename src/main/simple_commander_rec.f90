! concrete commander: 3D reconstruction routines
module simple_commander_rec
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
implicit none

public :: reconstruct3D_commander_distr
public :: reconstruct3D_commander
public :: volassemble_commander
public :: random_rec_commander_distr
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: reconstruct3D_commander_distr
  contains
    procedure :: execute      => exec_reconstruct3D_distr
end type reconstruct3D_commander_distr
type, extends(commander_base) :: reconstruct3D_commander
  contains
    procedure :: execute      => exec_reconstruct3D
end type reconstruct3D_commander
type, extends(commander_base) :: volassemble_commander
  contains
    procedure :: execute      => exec_volassemble
end type volassemble_commander
type, extends(commander_base) :: random_rec_commander_distr
  contains
    procedure :: execute      => exec_random_rec
end type random_rec_commander_distr

contains

    subroutine exec_reconstruct3D_distr( self, cline )
        class(reconstruct3D_commander_distr), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        character(len=LONGSTRLEN), allocatable :: list(:)
        character(len=:),          allocatable :: target_name
        character(len=STDLEN),     allocatable :: state_assemble_finished(:)
        type(reconstruct3D_commander)          :: xrec3D_shmem
        character(len=LONGSTRLEN) :: refine_path
        character(len=STDLEN)     :: volassemble_output, str_state, fsc_file
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        integer          :: state, ipart, sz_list
        logical          :: fall_over
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        if( .not. cline%defined('nparts') )then
            call xrec3D_shmem%execute(cline)
            return
        endif
        call build%init_params_and_build_spproj(cline, params)
        ! sanity check
        fall_over = .false.
        select case(trim(params%oritype))
            case('ptcl3D')
                fall_over = build%spproj%get_nptcls() == 0
            case('cls3D')
                fall_over = build%spproj%os_out%get_noris() == 0
        case DEFAULT
            THROW_HARD('unsupported ORITYPE')
        end select
        if( fall_over ) THROW_HARD('No images found!')
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        ! splitting
        if( trim(params%oritype).eq.'ptcl3D' ) call build%spproj%split_stk(params%nparts, dir=PATH_PARENT)
        ! eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! particle weights
        select case(trim(params%ptclw))
            case('yes')
                ! weights are set at search time, so nothing to do here.
            case DEFAULT
                call build%spproj_field%calc_hard_weights(params%frac)
        end select
        ! to update eo flags and weights
        call build%spproj%write_segment_inside(params%oritype)
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        ! assemble volumes
        ! this is for parallel volassemble over states
        allocate(state_assemble_finished(params%nstates) )
        do state = 1, params%nstates
            state_assemble_finished(state) = 'VOLASSEMBLE_FINISHED_STATE'//int2str_pad(state,2)
        enddo
        cline_volassemble = cline
        call cline_volassemble%set('prg', 'volassemble')
        ! parallel assembly
        do state = 1,params%nstates
            str_state = int2str_pad(state,2)
            volassemble_output = 'RESOLUTION_STATE'//trim(str_state)
            call cline_volassemble%set( 'state', real(state) )
            if( params%nstates>1 )call cline_volassemble%set('part', real(state))
            call qenv%exec_simple_prg_in_queue_async(cline_volassemble,&
            'simple_script_state'//trim(str_state), trim(volassemble_output))
        end do
        call qsys_watcher(state_assemble_finished)
        ! updates project file only if called from another workflow
        if( params%mkdir.eq.'yes' )then
            do state = 1,params%nstates
                fsc_file      = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                call build%spproj%add_fsc2os_out(trim(fsc_file), state, params%box)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd, state, 'vol_cavg')
                else
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd, state, 'vol')
                endif
            enddo
            call build%spproj%write_segment_inside('out',params%projfile)
        endif
        ! termination
        call qsys_cleanup
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D_distr

    subroutine exec_reconstruct3D( self, cline )
        class(reconstruct3D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        integer          :: s
        character(len=:), allocatable :: fbody
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_rec_eo_tbox(params)
        if( .not. cline%defined('nparts') )then ! shared-memory implementation
            ! eo partitioning
            if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
            ! particle weights
            select case(trim(params%ptclw))
                case('yes')
                ! weights are set at search time, so nothing to do here.
                case DEFAULT
                    call build%spproj_field%calc_hard_weights(params%frac)
            end select
            ! to update eo flags and weights
            call build%spproj%write_segment_inside(params%oritype)
        endif
        do s=1,params%nstates
            if( build%spproj_field%get_pop(s, 'state') == 0 ) cycle ! empty state
            fbody = 'recvol_state'
            call build%eorecvol%eorec(build%spproj, build%spproj_field, build%pgrpsyms, s, fbody=fbody)
        end do
        call qsys_job_finished( 'simple_rec_master :: exec_eorec')
        write(logfhandle,'(a,1x,a)') "GENERATED VOLUMES: reconstruct3D*.ext"
        ! end gracefully
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor_eo, only: reconstructor_eo
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
        logical                       :: l_euclid_regularization
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
        call build%eorecvol%kill_exp         ! reduced meory usage
        allocate(res05s(params%nstates), res0143s(params%nstates))
        res0143s = 0.
        res05s   = 0.
        call eorecvol_read%new( build%spproj)
        call eorecvol_read%kill_exp ! reduced memory usage
        n = params%nstates*params%nparts
        l_euclid_regularization = (params%cc_objfun==OBJFUN_EUCLID) .or. params%l_needs_sigma
        if( params%l_nonuniform ) l_euclid_regularization = .false.
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
        do ss=1,params%nstates
            if( cline%defined('state') )then
                s     = 1             ! index in reconstruct3D
                state = params%state  ! actual state
            else
                s     = ss
                state = ss
            endif
            if( L_BENCH_GLOB )then
                rt_read       = 0.
                rt_sum_reduce = 0.
            endif
            call build%eorecvol%reset_all
            ! assemble volumes
            do part=1,params%nparts
                if( L_BENCH_GLOB ) t_read = tic()
                call eorecvol_read%read_eos(trim(VOL_FBODY)//int2str_pad(state,2)//'_part'//int2str_pad(part,params%numlen))
                ! sum the Fourier coefficients
                if( L_BENCH_GLOB )then
                    rt_read       = rt_read + toc(t_read)
                    t_sum_reduce  = tic()
                endif
                call build%eorecvol%sum_reduce(eorecvol_read)
                if( L_BENCH_GLOB ) rt_sum_reduce = rt_sum_reduce + toc(t_sum_reduce)
            end do
            ! correct for sampling density and estimate resolution
            allocate(recname, source=trim(VOL_FBODY)//int2str_pad(state,2))
            allocate(volname, source=recname//params%ext)
            eonames(1) = trim(recname)//'_even'//params%ext
            eonames(2) = trim(recname)//'_odd'//params%ext
            resmskname = params%mskfile
            if( l_euclid_regularization )then
                ! the sum is done after regularization
            else
                if( L_BENCH_GLOB ) t_sum_eos = tic()
                call build%eorecvol%sum_eos
                if( L_BENCH_GLOB ) rt_sum_eos = rt_sum_eos + toc(t_sum_eos)
            endif
            if( L_BENCH_GLOB )then
                t_sampl_dens_correct_eos = tic()
            endif
            call build%eorecvol%sampl_dens_correct_eos(state, eonames(1), eonames(2), find4eoavg)
            if( L_BENCH_GLOB )then
                rt_sampl_dens_correct_eos = rt_sampl_dens_correct_eos + toc(t_sampl_dens_correct_eos)
            endif
            if( l_euclid_regularization )then
                if( L_BENCH_GLOB ) t_sum_eos = tic()
                call build%eorecvol%sum_eos
                if( L_BENCH_GLOB ) rt_sum_eos = rt_sum_eos + toc(t_sum_eos)
            endif
            call build%eorecvol%get_res(res05s(s), res0143s(s))
            if( L_BENCH_GLOB ) t_sampl_dens_correct_sum = tic()
            call build%eorecvol%sampl_dens_correct_sum( build%vol )
            if( L_BENCH_GLOB ) rt_sampl_dens_correct_sum = rt_sampl_dens_correct_sum + toc(t_sampl_dens_correct_sum)
            call build%vol%write( volname, del_if_exists=.true. )
            call wait_for_closure( volname )
            ! need to put the sum back at lowres for the eo pairs
            if( L_BENCH_GLOB ) t_eoavg = tic()
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
            if( L_BENCH_GLOB ) rt_eoavg = rt_eoavg + toc(t_eoavg)
            deallocate(recname, volname)
            if( cline%defined('state') )exit
        end do
        ! set the resolution limit according to the worst resolved model
        res  = maxval(res0143s)
        params%lp = max(params%lpstop, res)
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
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'VOLASSEMBLE_BENCH.txt'
            call fopen(fnr, FILE=trim(benchfname), STATUS='REPLACE', action='WRITE')
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

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander_distr), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(reconstruct3D_commander_distr) :: xrec3D_distr
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes'   )
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%os_ptcl3D%rnd_oris
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call cline%set('mkdir', 'no') ! to avoid nested dirs
        call cline%set('prg',   'reconstruct3D')
        call xrec3D_distr%execute(cline)
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_random_rec

end module simple_commander_rec
