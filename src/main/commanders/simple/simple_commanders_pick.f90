!@descr: for picking, extraction, and making picking references
module simple_commanders_pick
use simple_commanders_api
implicit none

public :: commander_pick
public :: commander_extract
public :: commander_reextract
public :: commander_pick_extract
public :: commander_make_pickrefs
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_pick
  contains
    procedure :: execute      => exec_pick
end type commander_pick

type, extends(commander_base) :: commander_extract
  contains
    procedure :: execute      => exec_extract
end type commander_extract

type, extends(commander_base) :: commander_reextract
  contains
    procedure :: execute      => exec_reextract
end type commander_reextract

type, extends(commander_base) :: commander_pick_extract
  contains
    procedure :: execute      => exec_pick_extract
end type commander_pick_extract

type, extends(commander_base) :: commander_make_pickrefs
  contains
    procedure :: execute      => exec_make_pickrefs
end type commander_make_pickrefs

contains


    subroutine exec_pick( self, cline )
        use simple_core_module_api, only: simple_end
        use simple_pick_strategy,   only: pick_strategy, create_pick_strategy
        use simple_cmdline,         only: cmdline
        use simple_parameters,      only: parameters
        class(commander_pick), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
        class(pick_strategy), allocatable :: strategy
        type(parameters) :: params
        call cline%set('prg', 'pick')
        strategy = create_pick_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_pick

    subroutine exec_extract( self, cline )
        use simple_core_module_api,   only: simple_end
        use simple_extract_strategy,  only: extract_strategy, create_extract_strategy
        use simple_cmdline,           only: cmdline
        use simple_parameters,        only: parameters
        class(commander_extract), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        class(extract_strategy), allocatable :: strategy
        type(parameters) :: params
        ! Helps distributed job script generation if it relies on 'prg'
        call cline%set('prg', 'extract')
        strategy = create_extract_strategy(cline)
        call strategy%apply_defaults(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_extract

    subroutine exec_reextract( self, cline )
        use simple_core_module_api,      only: simple_end
        use simple_reextract_strategy,   only: reextract_strategy, create_reextract_strategy
        use simple_cmdline,              only: cmdline
        use simple_parameters,           only: parameters
        class(commander_reextract), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        class(reextract_strategy), allocatable :: strategy
        type(parameters) :: params
        call cline%set('prg', 'reextract')
        strategy = create_reextract_strategy(cline)
        call strategy%initialize(params, cline)
        call strategy%execute(params, cline)
        call strategy%finalize_run(params, cline)
        call strategy%cleanup(params, cline)
        call simple_end(strategy%end_message())
        if( allocated(strategy) ) deallocate(strategy)
    end subroutine exec_reextract

    ! Stream only application
    subroutine exec_pick_extract( self, cline )
        use simple_picker_iter, only: picker_iter
        class(commander_pick_extract), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        type(parameters)              :: params
        type(oris)                    :: os_mic
        type(ori)                     :: o_mic
        type(picker_iter)             :: piter
        type(commander_extract)       :: xextract
        type(cmdline)                 :: cline_extract
        type(sp_project)              :: spproj
        type(string) :: micname, output_dir_picker, fbody, output_dir_extract
        type(string) :: boxfile, thumb_den
        integer :: fromto(2), imic, ntot, state, nvalid, i, nptcls
        logical :: l_extract
        ! set oritype
        call cline%set('oritype', 'mic')
        ! parse parameters
        call params%new(cline)
        ! if( params%stream.ne.'yes' ) THROW_HARD('new streaming only application')
        l_extract   = trim(params%extract).eq.'yes'
        ! read in movies
        call spproj%read( params%projfile )
        if( spproj%get_nintgs() == 0 ) THROW_HARD('no micrograph to process!')
        params%smpd = spproj%os_mic%get(1,'smpd')
        call cline%set('smpd',params%smpd)
        ! output directories
        output_dir_picker  = DIR_PICKER
        if( l_extract ) output_dir_extract = DIR_EXTRACT
        if( cline%defined('dir') )then
            output_dir_picker  = filepath(params%dir,output_dir_picker)//'/'
            if( l_extract) output_dir_extract = filepath(params%dir,output_dir_extract)//'/'
        endif
        call simple_mkdir(output_dir_picker)
        if( l_extract ) call simple_mkdir(output_dir_extract)
        ! picker specs
        select case(trim(params%picker))
            case('new')
                if(cline%defined('pickrefs'))then
                else
                    if( .not.cline%defined('moldiam') )then
                        THROW_HARD('MOLDIAM required for picker=new')
                    endif
                endif
            case('segdiam')
                if( .not.cline%defined('moldiam_max') )then
                    THROW_WARN('MOLDIAM_MAX not set on command line, falling back on default value: '//int2str(int(params%moldiam_max))//' A')
                endif
            case DEFAULT
                THROW_HARD('Unsupported PICKER: '//trim(params%picker))
        end select
        ! command lines
        if( l_extract )then
            cline_extract = cline
            call cline_extract%set('dir', output_dir_extract)
            call cline_extract%set('pcontrast', params%pcontrast)
            if( cline%defined('box_extract') ) call cline_extract%set('box', params%box_extract)
            call cline%delete('box')
            call cline_extract%delete('box_extract')
        endif
        ! file name
        if( cline%defined('fbody') )then
            fbody = params%fbody
        else
            fbody = ''
        endif
        ! range
        fromto(:) = 1
        if( cline%defined('fromp') .and. cline%defined('top') )then
            fromto(:) = [params%fromp, params%top]
        endif
        ntot   = fromto(2) - fromto(1) + 1
        nvalid = 0
        ! main loop
        do imic = fromto(1),fromto(2)
            ! fetch movie orientation
            call spproj%os_mic%get_ori(imic, o_mic)
            ! sanity check
            state = 1
            if( o_mic%isthere('state') ) state = o_mic%get_state()
            if( state == 0 ) cycle
            if( .not.o_mic%isthere('intg')   )cycle
            call o_mic%getter('intg', micname)
            if( .not.file_exists(micname)) cycle
            ! picker
            params%lp = max(params%fny, params%lp_pick)
            call piter%iterate(params, cline, params%smpd, micname, output_dir_picker, boxfile, thumb_den, nptcls)
            call o_mic%set('nptcls', nptcls)
            if( nptcls > 0 )then
                call o_mic%set('boxfile', boxfile)
                call o_mic%set('thumb_den', thumb_den)
            else
                call o_mic%set_state(0)
            endif
            ! update project
            call spproj%os_mic%set_ori(imic, o_mic)
            nvalid = nvalid+1
        enddo
        ! extract particles
        if( l_extract )then
            call spproj%write_segment_inside(params%oritype, params%projfile)
            call xextract%execute(cline_extract)
            ! nothing to write, done by extract
        else
            if( ntot > 1 )then
                ! purging state=0 and nptcls=0 mics such that all mics (nmics>1)
                ! can be assumed valid
                call os_mic%new(nvalid, is_ptcl=.false.)
                i = 0
                do imic = fromto(1),fromto(2)
                    state  = spproj%os_mic%get_state(imic)
                    nptcls = spproj%os_mic%get_int(imic,'nptcls')
                    if( (state == 1) .and. (nptcls > 0) )then
                        i = i+1
                        call os_mic%transfer_ori(i, spproj%os_mic, imic)
                    endif
                enddo
                spproj%os_mic = os_mic
                call os_mic%kill
            endif
            call spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        ! end gracefully
        call qsys_job_finished(params, string('simple_commanders_pick :: exec_pick_extract'))
        call o_mic%kill
        call piter%kill
        call simple_end('**** SIMPLE_PICK_EXTRACT NORMAL STOP ****')
    end subroutine exec_pick_extract

    subroutine exec_make_pickrefs( self, cline )
        use simple_pick_strategy,   only: make_pickrefs_impl
        use simple_core_module_api, only: simple_end
        class(commander_make_pickrefs), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call cline%set('prg', 'make_pickrefs')
        call make_pickrefs_impl(cline)
        call simple_end('**** SIMPLE_MAKE_PICKREFS NORMAL STOP ****')
    end subroutine exec_make_pickrefs

end module simple_commanders_pick
