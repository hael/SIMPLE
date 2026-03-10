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

type, extends(commander_base) :: commander_rec3D_worker
  contains
    procedure :: execute      => exec_rec3D_distr_worker
end type commander_rec3D_worker

type, extends(commander_base) :: random_rec_commander
  contains
    procedure :: execute      => exec_random_rec
end type random_rec_commander

contains

    subroutine exec_rec3D( self, cline )
        use simple_rec3D_strategy, only: rec3D_strategy, create_rec3D_strategy
        use simple_parameters,     only: parameters
        use simple_builder,        only: builder
        class(commander_rec3D), intent(inout) :: self
        class(cmdline),         intent(inout) :: cline
        class(rec3D_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        ! Commander-level defaults (apply to both modes)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)     ! to assure that shifts are being used
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        ! Select and run strategy
        strategy = create_rec3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        ! End gracefully (single unified termination)
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_rec3D

    subroutine exec_rec3D_distr_worker( self, cline )
        class(commander_rec3D_worker), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        type(string)         :: fname
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
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
    end subroutine exec_rec3D_distr_worker

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(parameters)      :: params
        type(builder)         :: build
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
