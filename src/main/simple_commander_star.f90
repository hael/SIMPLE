! concrete commander: Support operations on projects for STAR-formatted files
module simple_commander_star
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_star,           only: star_project
use simple_oris,           only: oris
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters
implicit none

public :: exportstar_project_commander
public :: importstar_project_commander
public :: print_star_project_info_commander
private

type, extends(commander_base) :: exportstar_project_commander
  contains
    procedure :: execute      => exec_exportstar_project
end type
type, extends(commander_base) :: importstar_project_commander
  contains
    procedure :: execute      => exec_importstar_project
end type
type, extends(commander_base) ::  print_star_project_info_commander
  contains
    procedure :: execute      => exec_print_star_project_info
end type
contains

    !> convert text (.txt) oris doc to binary (.simple)
    subroutine exec_exportstar_project( self, cline )
        class(exportstar_project_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(oris)       :: os
        type(sp_project) :: spproj
        integer          :: noris, ppos
        call params%new(cline)
        if( params%starfile == 'NONE')then
            ppos = scan(trim(params%projfile),".", BACK= .true.)
            if ( ppos <= 0 )then
                call simple_stop( " exportstar_project must be within a valid SIMPLE project. Projfile failed.")
            end if
            params%starfile = params%projfile(1:ppos)//'star'
        endif

        noris = nlines(params%oritab)
        call os%new(noris)
        call os%read(params%oritab)
        if( file_exists(params%projfile) )then
            call spproj%read(params%projfile)
        endif
        call spproj%set_sp_oris(params%oritype, os)
        ! Convert Simple's single-particle project file to STAR
        call spproj%write(params%projfile)


        call spproj%kill
        call simple_end('**** EXPORTSTAR_PROJECT NORMAL STOP ****')

    end subroutine exec_exportstar_project

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_importstar_project( self, cline )
        class(importstar_project_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline)
        if( params%starfile == 'NONE')&
            call simple_stop( " importstar_project must have input file, starfile=<filename>")
        if(.not. file_exists(params%starfile))&
            call simple_stop( " importstar_project must have a valid input file, starfile=<filename>")

        ! Import STAR filename

        ! Copy to SIMPLE's project format

        call spproj%read_segment(params%oritype, params%projfile)
        call spproj%write_segment2txt(params%oritype, params%outfile)
        call spproj%kill
        call simple_end('**** IMPORTSTAR_PROJECT NORMAL STOP ****')
    end subroutine exec_importstar_project


    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_print_star_project_info( self, cline )
        class(print_star_project_info_commander), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        call params%new(cline)
        call spproj%read(params%projfile)
        call spproj%print_info
        call spproj%kill
        call simple_end('**** PRINT_PROJECT_STAR_INFO NORMAL STOP ****')
    end subroutine exec_print_star_project_info



end module simple_commander_star
