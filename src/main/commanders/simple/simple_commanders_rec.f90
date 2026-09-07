!@descr: 3D reconstruction and associated things
module simple_commanders_rec
use simple_commanders_api
use simple_matcher_2Dprep
use simple_matcher_3Drec, only: calc_3Drec, calc_projdir3Drec
use simple_sigma2_files, only: load_sigma2_groups
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
        type(string)     :: rec_backend
        ! Commander-level defaults (apply to both modes)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir', 'yes')
        if( .not. cline%defined('trs')     ) call cline%set('trs', 5.)     ! to assure that shifts are being used
        if( .not. cline%defined('rec_backend') ) call cline%set('rec_backend', 'gridding')
        rec_backend = cline%get_carg('rec_backend')
        call cline%set('oritype', 'ptcl3D')
        call cline%delete('refine')
        ! Select and run strategy
        strategy = create_rec3D_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params, build, cline)
        call rec_backend%kill
        ! End gracefully (single unified termination)
        call simple_end('**** SIMPLE_RECONSTRUCT3D NORMAL STOP ****', print_simple=.false.)
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_rec3D

    subroutine exec_rec3D_distr_worker( self, cline )
        use simple_rec3D_pcg_strategy, only: execute_rec3D_pcg_worker
        class(commander_rec3D_worker), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)     :: params
        type(builder)        :: build
        integer, allocatable :: pinds(:)
        integer              :: nptcls2update
        logical              :: l_sigma_loaded
        call build%init_params_and_build_general_tbox(cline, params)
        call build%build_strategy3D_tbox(params)
        if( params%l_update_frac .and. build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls2update, pinds)
        else
            ! we sample all state > 0 and updatecnt > 0
            call build%spproj_field%sample4rec([params%fromp,params%top], nptcls2update, pinds)
        endif
        if( trim(params%rec_backend) == 'pcg' )then
            call execute_rec3D_pcg_worker(params, build, cline, pinds)
        else
            if( params%cc_objfun == OBJFUN_EUCLID )then
                call load_sigma2_groups(params, build%pftc, build%esig, build%spproj, build%spproj_field, &
                    &cline, l_sigma_loaded)
                if( .not. l_sigma_loaded ) THROW_HARD('gridding objfun=euclid requires sigma2 files')
            endif
            if( trim(params%projrec) == 'yes' )then
                call calc_projdir3Drec(params, build, cline, nptcls2update, pinds)
            else
                call calc_3Drec(params, build, cline, nptcls2update, pinds)
            endif
        endif
        ! cleanup
        call build%esig%kill
        call build%kill_strategy3D_tbox
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_rec :: exec_rec3D'))
    end subroutine exec_rec3D_distr_worker

    subroutine exec_random_rec( self, cline )
        class(random_rec_commander), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        type(commander_rec3D) :: xrec3D
        type(parameters)      :: params
        type(builder)         :: build
        call cline%set('oritype', 'ptcl3D')
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'yes'   )
        call build%init_params_and_build_spproj(cline, params)
        call build%spproj%os_ptcl3D%rnd_oris
        call build%spproj%write_segment_inside('ptcl3D', params%projfile)
        call cline%set('mkdir', 'no') ! to avoid nested dirs
        call cline%set('prg',   'rec3D')
        call xrec3D%execute(cline)
        call build%spproj_field%kill
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_RANDOM_REC NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_random_rec

end module simple_commanders_rec
