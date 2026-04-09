!@descr: reprojection commanders
module simple_commanders_reproject
use simple_commanders_api
use simple_simple_volinterp, only: reproject, rotvol
use simple_projector,        only: projector
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_reproject
 contains
   procedure :: execute      => exec_reproject
end type commander_reproject

contains
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

end module simple_commanders_reproject
