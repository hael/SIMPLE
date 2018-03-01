! concrete commander: operations on projects (sp_project) and associated files
module simple_commander_project
#include "simple_lib.f08"
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_params,         only: params

use simple_sp_project,     only: sp_project
implicit none

public :: project2txt_commander
public :: txt2project_commander
public :: print_project_info_commander
public :: manage_project_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: project2txt_commander
  contains
    procedure :: execute      => exec_project2txt
end type project2txt_commander
type, extends(commander_base) :: txt2project_commander
  contains
    procedure :: execute      => exec_txt2project
end type txt2project_commander
type, extends(commander_base) :: print_project_info_commander
  contains
    procedure :: execute      => exec_print_project_info
end type print_project_info_commander
type, extends(commander_base) :: manage_project_commander
  contains
    procedure :: execute      => exec_manage_project
end type manage_project_commander

contains

    !> convert text (.txt) oris doc to binary (.simple)
    subroutine exec_txt2project( self, cline )
        use simple_oris, only: oris
        class(txt2project_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)     :: p
        type(oris)       :: os
        type(sp_project) :: sp_proj
        integer          :: noris
        p     = params(cline)
        noris = nlines(p%oritab)
        call os%new_clean(noris)
        call os%read(p%oritab)
        if( file_exists(p%projfile) )then
            call sp_proj%read(p%projfile)
        endif
        call sp_proj%set_sp_oris(p%oritype, os)
        call sp_proj%write(p%projfile)
        call sp_proj%kill
        call simple_end('**** TXT2PROJECT NORMAL STOP ****')
    end subroutine exec_txt2project

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_project2txt( self, cline )
        class(project2txt_commander), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(params)     :: p
        type(sp_project) :: sp_proj
        p = params(cline)
        call sp_proj%read(p%projfile)
        call sp_proj%write_segment(p%oritype, p%outfile)
        call sp_proj%kill
        call simple_end('**** PROJECT2TXT NORMAL STOP ****')
    end subroutine exec_project2txt

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_print_project_info( self, cline )
        class(print_project_info_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(params)     :: p
        type(sp_project) :: sp_proj
        p = params(cline)
        call sp_proj%read(p%projfile)
        call sp_proj%print_info
        call sp_proj%kill
        call simple_end('**** PRINT_PROJECT_INFO NORMAL STOP ****')
    end subroutine exec_print_project_info

    !> for managing the projinfo and compenv segments of project
    subroutine exec_manage_project( self, cline )
        class(manage_project_commander), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(sp_project)      :: sp_proj
        character(len=STDLEN) :: projfile
        if( cline%defined('projfile') )then
            projfile = cline%get_carg('projfile')
            if( .not. file_exists(projfile) )then
                write(*,*) 'The inputted project file (projfile): ', trim(projfile)
                stop 'does not exist in cwd; commander_project :: exec_manage_project'
                call sp_proj%read(projfile)
            endif
        endif
        call sp_proj%update_projinfo( cline )
        call sp_proj%update_compenv( cline )
        call sp_proj%write
        ! debug
        call sp_proj%projinfo%write('projinfo.txt')
        call sp_proj%compenv%write('compenv.txt')
        call simple_end('**** PRINT_PROJECT_INFO NORMAL STOP ****')
    end subroutine exec_manage_project

end module simple_commander_project
