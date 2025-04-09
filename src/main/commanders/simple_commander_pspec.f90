module simple_commander_pspec
include 'simple_lib.f08'
use simple_commander_base,      only: commander_base
use simple_image,               only: image
use simple_pspecs,              only: pspecs
use simple_strategy2D3D_common, only: read_imgbatch
use simple_cmdline,             only: cmdline
use simple_parameters,          only: parameters, params_glob
use simple_sp_project,          only: sp_project
use simple_builder,             only: builder, build_glob
use simple_stack_io,            only: stack_io
implicit none

public :: analyze_pspecs_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: analyze_pspecs_commander
    contains
    procedure :: execute => exec_analyze_pspecs
end type analyze_pspecs_commander

contains

    subroutine exec_analyze_pspecs( self, cline )
        class(analyze_pspecs_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters) :: params
        type(builder)    :: build
        type(pspecs)     :: pows
        integer          :: nptcls, iptcl, ind_in_stk
        type(stack_io)   :: stkio_r
        character(len=:), allocatable :: stkname
        if( .not. cline%defined('hp')      ) call cline%set('hp',       20.)
        if( .not. cline%defined('lp')      ) call cline%set('lp',        6.)
        if( .not. cline%defined('mkdir')   ) call cline%set('mkdir',  'yes')
        call build%init_params_and_build_general_tbox(cline, params)
        call build%spproj%update_projinfo(cline)
        call build%spproj%write_segment_inside('projinfo')
        nptcls = build%spproj%get_nptcls()
        call pows%new(nptcls, build%img, params%hp, params%lp )
        do iptcl = 1, nptcls
            ! read image
            call build_glob%spproj%get_stkname_and_ind(params%oritype, iptcl, stkname, ind_in_stk)
            if( .not. stkio_r%stk_is_open() )then
                call stkio_r%open(stkname, params_glob%smpd, 'read')
            else if( .not. stkio_r%same_stk(stkname, [params%box,params_glob%box,1]) )then
                call stkio_r%close
                call stkio_r%open(stkname, params%smpd, 'read')
            endif
            call stkio_r%read(ind_in_stk, build%img)
            ! update pspec object
            call pows%set_pspec(iptcl, build%img, params%msk)
        end do
        call stkio_r%close
        call simple_end('**** SIMPLE_ANALYZE_PSPECS NORMAL STOP ****')
    end subroutine exec_analyze_pspecs

end module simple_commander_pspec