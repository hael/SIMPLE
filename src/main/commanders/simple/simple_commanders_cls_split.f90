!@descr: class split commander
module simple_commanders_cls_split
use simple_commanders_api
implicit none

type, extends(commander_base) :: commander_cls_split
  contains
    procedure :: execute => exec_cls_split
end type commander_cls_split

type, extends(commander_base) :: commander_diffmap_denoise_project
  contains
    procedure :: execute => exec_diffmap_denoise_project
end type commander_diffmap_denoise_project

contains

    subroutine exec_cls_split( self, cline )
        use simple_cls_split_strategy
        class(commander_cls_split), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(cls_split_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_cls_split_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_CLS_SPLIT NORMAL STOP ****')
    end subroutine exec_cls_split

    subroutine exec_diffmap_denoise_project( self, cline )
        use simple_diffmap_denoise_project_strategy
        class(commander_diffmap_denoise_project), intent(inout) :: self
        class(cmdline),                            intent(inout) :: cline
        class(diffmap_denoise_project_strategy), allocatable :: strategy
        type(parameters) :: params
        type(builder)    :: build
        strategy = create_diffmap_denoise_project_strategy(cline)
        call strategy%initialize(params, build, cline)
        call strategy%execute(params, build, cline)
        call strategy%finalize_run(params, build, cline)
        call strategy%cleanup(params)
        if( allocated(strategy) ) deallocate(strategy)
        call build%kill_general_tbox
        call simple_end('**** SIMPLE_DIFFMAP_DENOISE_PROJECT NORMAL STOP ****')
    end subroutine exec_diffmap_denoise_project

end module simple_commanders_cls_split
