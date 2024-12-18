! concrete commander: 3D reconstruction routines
module simple_commander_rec
include 'simple_lib.f08'
use simple_builder,        only: builder
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_parameters,     only: parameters
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs
use simple_strategy2D3D_common
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
        type(reconstruct3D_commander)          :: xrec3D_shmem
        character(len=STDLEN)       :: str_state, fsc_file, nthr_here_str
        type(volassemble_commander) :: xvolassemble
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        integer          :: state, envlen, io_stat, nthr_here
        logical          :: fall_over
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.) ! to assure that shifts are being used
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        if( .not. cline%defined('nparts') )then
            call xrec3D_shmem%execute(cline)
            return
        endif
         ! deal with # threads for the master process
        call get_environment_variable('SLURM_CPUS_PER_TASK', nthr_here_str, envlen)
        if( envlen > 0 )then
            call str2int(trim(nthr_here_str), io_stat, nthr_here)
        else
            !$ nthr_here = omp_get_max_threads()
            nthr_here = min(NTHR_SHMEM_MAX,nthr_here)
        end if
        write(logfhandle,'(A,I6)')'>>> # SHARED-MEMORY THREADS USED BY RECONSTRUCT3D MASTER PROCESS: ', nthr_here
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
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
        if( .not.(trim(params%ptclw).eq.'yes') ) call build%spproj_field%calc_hard_weights(params%frac)
        ! to update eo flags and weights
        call build%spproj%write_segment_inside(params%oritype)
        ! setup the environment for distributed execution
        call qenv%new(params%nparts)
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR)
        ! assemble volumes
        cline_volassemble = cline
        call cline_volassemble%set('prg', 'volassemble')
        call cline_volassemble%set('nthr',    nthr_here)
        call xvolassemble%execute_safe(cline_volassemble)
        ! updates project file only if mkdir is set to yes
        if( params%mkdir.eq.'yes' )then
            do state = 1,params%nstates
                str_state = int2str_pad(state,2)
                fsc_file  = FSC_FBODY//trim(str_state)//trim(BIN_EXT)
                call build%spproj%add_fsc2os_out(trim(fsc_file), state, params%box_crop)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd_crop, state, 'vol_cavg')
                else
                    call build%spproj%add_vol2os_out(trim(VOL_FBODY)//trim(str_state)//params%ext, params%smpd_crop, state, 'vol')
                endif
            enddo
            call build%spproj%write_segment_inside('out',params%projfile)
        endif
        ! termination
        call qsys_cleanup
        call build%spproj_field%kill
        call qenv%kill
        call cline_volassemble%kill
        call job_descr%kill
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D_distr

    subroutine exec_reconstruct3D( self, cline )
        use simple_euclid_sigma2, only: euclid_sigma2, eucl_sigma2_glob
        class(reconstruct3D_commander), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)                :: params
        type(builder)                   :: build
        type(euclid_sigma2),     target :: eucl_sigma
        character(len=:),   allocatable :: fname
        integer,            allocatable :: pinds(:), updatecnts(:)
        logical,            allocatable :: ptcl_mask(:)
        integer :: nptcls2update
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( .not. cline%defined('nparts') )then ! shared-memory implementation
            ! eo partitioning
            if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
            ! particle weights
            if( .not.(trim(params%ptclw).eq.'yes') ) call build%spproj_field%calc_hard_weights(params%frac)
            ! to update eo flags and weights
            call build%spproj%write_segment_inside(params%oritype)
        endif
        allocate(ptcl_mask(params%fromp:params%top)) 
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds, ptcl_mask)
        else
            ! we sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds, ptcl_mask)
        endif
        if( params%l_needs_sigma )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call eucl_sigma%new(fname, params%box)
            call eucl_sigma%read_groups(build%spproj_field, ptcl_mask)
        end if
        if( trim(params%projrec).eq.'yes' )then
            ! making sure the projection directions assignment
            ! refers to current reference space
            call build%spproj_field%set_projs(build%eulspace)
            call calc_projdir3Drec( cline, nptcls2update, pinds )
        else
            call calc_3Drec( cline, nptcls2update, pinds )
        endif
        ! cleanup
        call eucl_sigma%kill
        call qsys_job_finished('simple_commander_rec :: exec_reconstruct3D')
        ! end gracefully
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reconstruct3D

    subroutine exec_volassemble( self, cline )
        use simple_reconstructor_eo, only: reconstructor_eo
        use simple_image,            only: image
        class(volassemble_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)              :: params
        type(builder)                 :: build
        type(reconstructor_eo)        :: eorecvol_read
        type(image)                   :: vol_prev_even, vol_prev_odd
        character(len=:), allocatable :: finished_fname, recname, volname, volname_prev, fsc_txt_file
        character(len=:), allocatable :: volname_prev_even, volname_prev_odd, str_state, str_iter
        character(len=LONGSTRLEN)     :: eonames(2), benchfname
        real, allocatable             :: fsc(:), res05s(:), res0143s(:)
        real                          :: weight_prev, update_frac_trail_rec
        integer                       :: part, state, find4eoavg, fnr
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
        call eorecvol_read%new(build%spproj, expand=.false.)
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
            call build%spproj%read_segment(params%oritype, params%projfile)
            update_frac_trail_rec = build%spproj%os_ptcl3D%get_update_frac()
        endif
        ! assemble volumes
        do state=1,params%nstates
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
            allocate(recname, source=VOL_FBODY//int2str_pad(state,2))
            allocate(volname, source=recname//params%ext)
            eonames(1) = trim(recname)//'_even'//params%ext
            eonames(2) = trim(recname)//'_odd'//params%ext
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
                if( .not. file_exists(volname_prev_even) ) THROW_HARD('File: '//trim(volname_prev_even)//' does not exist!')
                if( .not. file_exists(volname_prev_odd)  ) THROW_HARD('File: '//trim(volname_prev_odd)//' does not exist!')
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
                fsc_txt_file = 'RESOLUTION_STATE'//str_state//'_ITER'//str_iter
            else
                fsc_txt_file = 'RESOLUTION_STATE'//str_state
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
            call build%vol%write( volname, del_if_exists=.true. )
            call wait_for_closure( volname )
            ! need to put the sum back at lowres for the eo pairs
            if( L_BENCH_GLOB ) t_eoavg = tic()
            call build%vol%fft
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
            if( params%l_trail_rec .and. update_frac_trail_rec < 0.99 )then
                call build%vol%ifft
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
            deallocate(recname, volname)
        end do
        call eorecvol_read%kill
        ! end gracefully
        call simple_end('**** SIMPLE_VOLASSEMBLE NORMAL STOP ****', print_simple=.false.)
        ! indicate completion (when run in a qsys env)
        allocate( finished_fname, source='VOLASSEMBLE_FINISHED' )
        call simple_touch( finished_fname , errmsg='In: commander_rec::volassemble')
        if( L_BENCH_GLOB )then
            rt_tot     = toc(t_tot)
            benchfname = 'VOLASSEMBLE_BENCH_ITER'//int2str_pad(params%which_iter,3)//'.txt'
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
