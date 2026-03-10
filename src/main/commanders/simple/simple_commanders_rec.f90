!@descr: 3D reconstruction and associated things
module simple_commanders_rec
use simple_commanders_api
use simple_strategy2D3D_common
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_rec3D
  contains
    procedure :: execute => exec_rec3D
end type commander_rec3D

type, extends(commander_base) :: commander_rec3D_distr
  contains
    procedure :: execute => exec_rec3D_distr
end type commander_rec3D_distr

! type, extends(commander_base) :: commander_rec3D_distr_worker
!   contains
!     procedure :: execute      => exec_rec3D_distr_worker
! end type commander_rec3D_distr_worker

type, extends(commander_base) :: random_rec_commander
  contains
    procedure :: execute      => exec_random_rec
end type random_rec_commander

contains

    subroutine exec_rec3D_distr( self, cline )
        use simple_commanders_rec_distr, only: commander_volassemble
        class(commander_rec3D_distr), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(commander_rec3D) :: xrec3D_shmem
        type(commander_volassemble)   :: xvolassemble
        type(string)     :: str_state, fsc_file
        type(parameters) :: params
        type(builder)    :: build
        type(qsys_env)   :: qenv
        type(cmdline)    :: cline_volassemble
        type(chash)      :: job_descr
        integer          :: state, nthr_here
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
        call set_master_num_threads(nthr_here, string('rec3D'))
        ! parse parameters and project
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
        if( trim(params%oritype).eq.'ptcl3D' ) call build%spproj%split_stk(params%nparts, dir=string(PATH_PARENT))
        ! eo partitioning
        if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
        ! weights
        if( trim(params%ptclw).eq.'yes' )then
            ! not implemented
        else
            if( trim(params%cavgw).eq.'yes' )then
                ! class averages
                call build%spproj_field%calc_cavg_soft_weights(params%frac)
            else
                ! particles
                call build%spproj_field%calc_hard_weights(params%frac)
            endif
        endif
        ! to update eo flags and weights
        call build%spproj%write_segment_inside(params%oritype)
        ! setup the environment for distributed execution
        call qenv%new(params, params%nparts)
        call cline%gen_job_descr(job_descr)
        ! schedule
        call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        ! assemble volumes
        cline_volassemble = cline
        call cline_volassemble%set('prg', 'volassemble')
        call cline_volassemble%set('nthr',    nthr_here)
        call xvolassemble%execute(cline_volassemble)
        ! updates project file only if mkdir is set to yes
        if( params%mkdir.eq.'yes' )then
            do state = 1,params%nstates
                str_state = int2str_pad(state,2)
                fsc_file  = FSC_FBODY//str_state%to_char()//BIN_EXT
                call build%spproj%add_fsc2os_out(fsc_file, state, params%box_crop)
                if( trim(params%oritype).eq.'cls3D' )then
                    call build%spproj%add_vol2os_out(string(VOL_FBODY)//str_state//params%ext, params%smpd_crop, state, 'vol_cavg')
                else
                    call build%spproj%add_vol2os_out(string(VOL_FBODY)//str_state//params%ext, params%smpd_crop, state, 'vol')
                endif
            enddo
            call build%spproj%write_segment_inside('out',params%projfile)
        endif
        ! termination
        call qsys_cleanup(params)
        call build%spproj_field%kill
        call qenv%kill
        call cline_volassemble%kill
        call job_descr%kill
        call simple_end('**** SIMPLE_rec3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rec3D_distr

    subroutine exec_rec3D( self, cline )
        class(commander_rec3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(parameters)            :: params
        type(builder)               :: build
        type(string)                :: fname
        integer, allocatable        :: pinds(:)
        integer                     :: nptcls2update
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( .not. cline%defined('nparts') )then ! shared-memory implementation
            ! eo partitioning
            if( build%spproj_field%get_nevenodd() == 0 ) call build%spproj_field%partition_eo
            ! weights
            if( trim(params%ptclw).eq.'yes' )then
                ! not implemented
            else
                if( trim(params%cavgw).eq.'yes' )then
                    ! class averages
                    call build%spproj_field%calc_cavg_soft_weights(params%frac)
                else
                    ! particles
                    call build%spproj_field%calc_hard_weights(params%frac)
                endif
            endif
            ! to update eo flags and weights
            call build%spproj%write_segment_inside(params%oritype)
        endif
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            ! we sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        if( params%l_ml_reg )then
            fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
            call build%esig%new(params, build%pftc, fname, params%box)
            call build%esig%read_groups(build%pftc, build%spproj_field)
        end if
        call calc_3Drec( params, build, cline, nptcls2update, pinds )
        ! cleanup
        call build%esig%kill
        call qsys_job_finished(params, string('simple_commanders_rec :: exec_rec3D'))
        ! end gracefully
        call simple_end('**** SIMPLE_rec3D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_rec3D

!     subroutine exec_rec3D( self, cline )
!         use simple_rec3D_strategy, only: rec3D_strategy, create_rec3D_strategy
!         use simple_parameters,     only: parameters
!         use simple_builder,        only: builder
!         class(commander_rec3D), intent(inout) :: self
!         class(cmdline),         intent(inout) :: cline
!         class(rec3D_strategy), allocatable :: strategy
!         type(parameters) :: params
!         type(builder)    :: build
!         ! Commander-level defaults (apply to both modes)
!         if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
!         if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)     ! to assure that shifts are being used
!         if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
!         call cline%delete('refine')
!         ! Select and run strategy
!         strategy = create_rec3D_strategy(cline)
!         call strategy%initialize(params, build, cline)
!         call strategy%execute(params, build, cline)
!         call strategy%finalize_run(params, build, cline)
!         call strategy%cleanup(params, build, cline)
!         ! End gracefully (single unified termination)
!         call simple_end('**** SIMPLE_rec3D NORMAL STOP ****', print_simple=.false.)
!         if( allocated(strategy) ) deallocate(strategy)
!     end subroutine exec_rec3D

!    subroutine exec_rec3D_distr( self, cline )
!         class(commander_rec3D_distr), intent(inout) :: self
!         class(cmdline),               intent(inout) :: cline
!         type(commander_rec3D) :: xrec
!         call xrec%execute(cline)
!     end subroutine exec_rec3D_distr

!     subroutine exec_rec3D_distr_worker( self, cline )
!         class(commander_rec3D_distr_worker), intent(inout) :: self
!         class(cmdline),                      intent(inout) :: cline
!         type(parameters)     :: params
!         type(builder)        :: build
!         type(string)         :: fname
!         integer, allocatable :: pinds(:)
!         integer              :: nptcls2update
!         call build%init_params_and_build_general_tbox(cline, params)
!         call build%build_strategy3D_tbox(params)
!         if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
!             call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
!         else
!             ! we sample all state > 0 and updatecnt > 0
!             call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
!         endif
!         if( params%l_ml_reg )then
!             fname = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
!             call build%esig%new(params, build%pftc, fname, params%box)
!             call build%esig%read_groups(build%pftc, build%spproj_field)
!         end if
!         call calc_3Drec( params, build, cline, nptcls2update, pinds )
!         ! cleanup
!         call build%esig%kill
!         call qsys_job_finished(params, string('simple_commanders_rec :: exec_rec3D'))
!     end subroutine exec_rec3D_distr_worker

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(parameters) :: params
        type(builder)    :: build
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes'   )
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%os_ptcl3D%rnd_oris
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call cline%set('mkdir', 'no') ! to avoid nested dirs
        call cline%set('prg',   'rec3D')
        call xrec3D%execute(cline)
        call build%spproj_field%kill
        call simple_end('**** SIMPLE_RANDOM_REC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_random_rec

end module simple_commanders_rec
