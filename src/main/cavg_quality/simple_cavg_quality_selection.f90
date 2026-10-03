!@descr: scores a project's class averages with a quality model and writes the selected/rejected stacks
!==============================================================================
! MODULE: simple_cavg_quality_selection
!
! PURPOSE:
!   The two steps every class-average selection shares, whatever it then does
!   with the result:
!     - score_project_cavgs reads the class averages of a project and scores
!       them, with a quality model and relation parameters built from the
!       class-average box and mask diameter, or with the hard gates only;
!     - write_cavg_selection_stacks writes the selected and the rejected
!       class averages as two stacks.
!   The particle sieve (two tiers, rejection reasons, persistent compatibility
!   models), the model_cavgs_rejection commander (annotation, pruning, report
!   tables) and the stream's initial analysis (a fresh compatibility model)
!   keep their own policy and call these for the shared part.
!==============================================================================
module simple_cavg_quality_selection
use simple_defs,                  only: logfhandle
use simple_error,                 only: simple_exception
use simple_string,                only: string
use simple_fileio,                only: file_exists, del_file
use simple_image,                 only: image
use simple_cmdline,               only: cmdline
use simple_parameters,            only: parameters
use simple_sp_project,            only: sp_project
use simple_imgarr_utils,          only: read_cavgs_into_imgarr, dealloc_imgarr
use simple_cavg_quality_model,    only: cavg_quality_model
use simple_cavg_quality_types,    only: cavg_quality_result
use simple_cavg_quality_analysis, only: evaluate_cavg_quality, evaluate_cavg_quality_hard_reject
implicit none

public :: score_project_cavgs, write_cavg_selection_stacks, write_cavg_stack
private
#include "simple_local_flags.inc"

contains

    !> Reads the class averages of @p spproj into @p cavg_imgs and scores them into @p quality
    !! with @p model. With @p hard_gate_context only the hard gates of that context are applied.
    !! @p smpd is the class averages' pixel size for the relation parameters, which turn
    !! @p mskdiam into the mask radius (pixels) of the relational feature's signal statistics, as
    !! the model_cavgs_rejection commander (whose tables train the models) has it from the
    !! project. Without it they take the parameters default pixel size; that is harmless only
    !! with @p mskdiam = 0, where the radius falls back to the box default (the sieve).
    subroutine score_project_cavgs( spproj, model, mskdiam, cavg_imgs, quality, hard_gate_context, smpd )
        class(sp_project),          intent(inout) :: spproj
        type(cavg_quality_model),   intent(in)    :: model
        real,                       intent(in)    :: mskdiam
        type(image), allocatable,   intent(inout) :: cavg_imgs(:)
        type(cavg_quality_result),  intent(inout) :: quality
        character(len=*), optional, intent(in)    :: hard_gate_context
        real,             optional, intent(in)    :: smpd
        type(parameters) :: relation_params
        type(cmdline)    :: relation_cline
        if( allocated(cavg_imgs) ) call dealloc_imgarr(cavg_imgs)
        cavg_imgs = read_cavgs_into_imgarr(spproj)
        if( size(cavg_imgs) < 1 ) THROW_HARD('no class averages in the project; score_project_cavgs')
        if( present(hard_gate_context) )then
            call evaluate_cavg_quality_hard_reject(cavg_imgs, spproj%os_cls2D, mskdiam, quality, hard_gate_context)
            return
        endif
        call relation_cline%set('box',     cavg_imgs(1)%get_box())
        call relation_cline%set('oritype', 'cls2D')
        call relation_cline%set('ctf',     'no')
        call relation_cline%set('objfun',  'cc')
        call relation_cline%set('mskdiam', mskdiam)
        if( present(smpd) ) call relation_cline%set('smpd', smpd)
        call relation_params%new(relation_cline)
        call evaluate_cavg_quality(cavg_imgs, spproj%os_cls2D, mskdiam, quality, model, relation_params=relation_params)
        call relation_cline%kill
    end subroutine score_project_cavgs

    !> Writes the class averages with @p states > 0 to @p selected_fname and the others to
    !! @p rejected_fname, in class order, replacing existing files.
    subroutine write_cavg_selection_stacks( cavg_imgs, states, selected_fname, rejected_fname )
        class(image),  intent(inout) :: cavg_imgs(:)
        integer,       intent(in)    :: states(:)
        class(string), intent(in)    :: selected_fname, rejected_fname
        if( size(states) /= size(cavg_imgs) ) THROW_HARD('# states /= # class averages; write_cavg_selection_stacks')
        call write_cavg_stack(cavg_imgs, states > 0,         selected_fname)
        call write_cavg_stack(cavg_imgs, .not. (states > 0), rejected_fname)
    end subroutine write_cavg_selection_stacks

    !> Writes the class averages with @p mask set to @p fname, in class order, replacing an existing file.
    subroutine write_cavg_stack( cavg_imgs, mask, fname )
        class(image),  intent(inout) :: cavg_imgs(:)
        logical,       intent(in)    :: mask(:)
        class(string), intent(in)    :: fname
        integer :: icls, istk
        if( size(mask) /= size(cavg_imgs) ) THROW_HARD('mask size /= # class averages; write_cavg_stack')
        if( file_exists(fname) ) call del_file(fname)
        istk = 0
        do icls = 1,size(cavg_imgs)
            if( .not. mask(icls) ) cycle
            istk = istk + 1
            call cavg_imgs(icls)%write(fname, istk)
        enddo
        write(logfhandle,'(A,A,A,I6)') '>>> WROTE ', fname%to_char(), ' #CAVGS: ', istk
    end subroutine write_cavg_stack

end module simple_cavg_quality_selection
