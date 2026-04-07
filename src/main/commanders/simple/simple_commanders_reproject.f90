!@descr: reprojection commanders
module simple_commanders_reproject
use simple_commanders_api
use simple_simple_volinterp, only: reproject, rotvol
use simple_projector,        only: projector
implicit none
public :: run_reproj_polar_with_strategy
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_reproject
 contains
   procedure :: execute      => exec_reproject
end type commander_reproject

type, extends(commander_base) :: commander_reproj_polar
    contains
        procedure :: execute => exec_reproj_polar
end type commander_reproj_polar

contains

    subroutine run_reproj_polar_with_strategy( params, build, cline, runtime_ready, use_distr )
        use simple_reproj_polar_strategy, only: reproj_polar_strategy, reproject_distr_strategy, reproject_inmem_strategy, create_reproj_polar_strategy
        use simple_parameters,            only: parameters
        use simple_builder,               only: builder
        class(parameters), intent(inout) :: params
        class(builder),    intent(inout) :: build
        class(cmdline),    intent(inout) :: cline
        logical,           intent(in)    :: runtime_ready
        logical, optional, intent(in)    :: use_distr
        class(reproj_polar_strategy), allocatable :: strategy
        type(cmdline) :: cline_local
        if( runtime_ready )then
            cline_local = cline
            call cline_local%delete('part')
            call cline_local%delete('fromp')
            call cline_local%delete('top')
            strategy = allocate_strategy(cline_local, use_distr)
            call strategy%execute(params, build, cline_local)
            call strategy%cleanup(params, build, cline_local)
            call cline_local%kill
        else
            if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',   'no')
            if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl3D')
            if( .not. cline%defined('fromp')   ) call cline%set('fromp', 1)
            if( .not. cline%defined('top') .and. cline%defined('nspace') )then
                call cline%set('top', cline%get_iarg('nspace'))
            else if( .not. cline%defined('nspace') .and. cline%defined('top') )then
                call cline%set('nspace', cline%get_iarg('top'))
            else if( .not. cline%defined('nspace') .and. .not. cline%defined('top') )then
                THROW_HARD('Need either nspace or top for reproj_polar')
            endif
            strategy = allocate_strategy(cline, use_distr)
            call strategy%initialize(params, build, cline)
            call strategy%execute(params, build, cline)
            call strategy%cleanup(params, build, cline)
        endif
        if( allocated(strategy) ) deallocate(strategy)

        contains

            function allocate_strategy( cline_here, use_distr_here ) result(strategy_out)
                class(cmdline),        intent(in) :: cline_here
                logical,     optional, intent(in) :: use_distr_here
                class(reproj_polar_strategy), allocatable :: strategy_out
                if( present(use_distr_here) )then
                    if( use_distr_here )then
                        allocate(reproject_distr_strategy :: strategy_out)
                    else
                        allocate(reproject_inmem_strategy :: strategy_out)
                    endif
                else
                    strategy_out = create_reproj_polar_strategy(cline_here)
                endif
            end function allocate_strategy
            
    end subroutine run_reproj_polar_with_strategy

    !> exec_project generate projections from volume
    subroutine exec_reproject( self, cline )
        use simple_imgarr_utils, only: dealloc_imgarr
        class(commander_reproject), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(parameters)         :: params
        type(builder)            :: build
        type(image), allocatable :: imgs(:)
        integer,     allocatable :: states(:), tmp(:)
        integer                  :: i, s
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('outstk') ) call cline%set('outstk', 'reprojs'//STK_EXT)
        if( .not. cline%defined('oritab') )then
            if( .not. cline%defined('nspace') ) THROW_HARD('need nspace (for number of projections)!')
        endif
        call params%new(cline)
        if( cline%defined('oritab') )then
            params%nptcls = binread_nlines(params%oritab, params%spproj_iseg)
            call build%build_spproj(params, cline)
            call build%build_general_tbox(params, cline)
            params%nspace = build%spproj_field%get_noris()
        else
            params%nptcls = params%nspace
            call build%build_spproj(params, cline)
            call build%build_general_tbox(params, cline)
            call build%pgrpsyms%build_refspiral(build%spproj_field)
            call build%spproj_field%set_all2single('state',1)
        endif
        ! generate projections
        if( file_exists(params%outstk) ) call del_file(params%outstk)
        states = nint(build%spproj_field%get_all('state'))
        tmp    = states
        do s = 1,params%nstates
            if( cline%defined('state') )then
                if( s /= params%state ) cycle
            endif
            if( .not.cline%defined('vol'//int2str(s)) ) cycle
            if( any(states==s) )then
                ! read and mask
                call build%vol%read(params%vols(s))
                if(cline%defined('mskdiam')) call build%vol%mask3D_soft(params%msk,backgr=0.)
                ! project
                where( states == s )
                    tmp = 1
                elsewhere
                    tmp = 0
                end where
                call build%spproj_field%set_all('state', tmp)
                imgs = reproject(build%vol, build%spproj_field)
                ! write
                do i = 1,params%nspace
                    if( states(i) /= s ) cycle
                    if( params%neg .eq. 'yes' ) call imgs(i)%neg()
                    call imgs(i)%write(params%outstk,i)
                enddo
            endif
        enddo
        call build%spproj_field%set_all('state', states) ! restore
        call build%spproj_field%write(string('reproject_oris'//TXT_EXT), [1,params%nptcls])
        if( cline%defined('state') ) where( states /= params%state ) states = 0
        ! zero state=0
        call imgs(1)%zero
        do i = 1,params%nspace
            if( states(i) == 0 ) call imgs(1)%write(params%outstk,i)
        enddo
        call update_stack_nimgs(params%outstk, params%nspace)
        ! cleanup
        call build%kill_general_tbox
        call dealloc_imgarr(imgs)
        deallocate(tmp,states)
        call simple_end('**** SIMPLE_REPROJECT NORMAL STOP ****')
    end subroutine exec_reproject

    subroutine exec_reproj_polar( self, cline )
        use simple_parameters, only: parameters
        use simple_builder,    only: builder
        class(commander_reproj_polar), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        call run_reproj_polar_with_strategy(params, build, cline, .false.)
        call simple_end('**** SIMPLE_POLAR_REPROJ NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_reproj_polar

end module simple_commanders_reproject
