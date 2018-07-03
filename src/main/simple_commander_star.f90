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
        type(parameters)  :: params
        type(oris)        :: os
        type(sp_project)  :: spproj
        type(star_project):: starproj
        character(len=:),allocatable :: this_starfile
        integer           :: noris, ppos,iostat
        call params%new(cline)
        if( params%starfile == 'NONE')then
            ppos = scan(trim(params%projfile),".", BACK= .true.)
            if ( ppos <= 0 )then
                call simple_stop( " exportstar_project must be within a valid SIMPLE project. Projfile failed.")
            end if
            allocate(this_starfile, source=params%projfile(1:ppos)//'star')
        else
            allocate(this_starfile, source=trim(params%starfile))
        end if
        ! create backup in case the starfile exists
        if(file_exists(this_starfile))then 
            iostat= simple_rename(trim(this_starfile),&
                trim(this_starfile)//"-bkup")
        end if

        if( params%export_type  == 'NONE')then
           write (*,*) " exportstar_project must have a valid export_type."
           write (*,*) " Accepted values are micrographs|mcmicrographs|ctf_estimation|select|extract|class2d|init3dmodel|refine3d|post or all."
           write (*,*) "    micrographs:    for non-dose weighted micrographs"
           write (*,*) "    mcmicrographs:  for motion-corrected, dose weighted micrographs"
           write (*,*) "    ctf_estimation: for ctf estimation and micrographs"

           write (*,*) "    class2d:   for 2D class averages"
           write (*,*) "    select:    for selected 2D class averages"

           write (*,*) "    init3dmodel: for initial particles"
           write (*,*) "    extract:   for extract dose-weighted particles"
           write (*,*) "    refine3d:  for refined 3D particles"
           write (*,*) "    post:      for post-processing"

 call simple_stop( " exportstar_project failed due to incorrect export_type arg")
end if
        !! Prepare project
        call starproj%prepare(spproj,  this_starfile)
     
        ! Use export type to define export output variables
        select case(params%export_type)
        case('micrographs')
           call  starproj%export_micrographs(spproj,  params%starfile)
        case('mcmicrographs')
           call  starproj%export_motion_corrected_micrographs(spproj,  params%starfile)
        case('ctf_estimation')
           call starproj%export_ctf_estimation(spproj,  params%starfile)
        case('select')
            call starproj%export_class2D_select(spproj,  params%starfile)
        case('extract')
           call starproj%export_extract_doseweightedptcls(spproj,  params%starfile)
        case('class2d')
           call starproj%export_class2D(spproj,  params%starfile)
        case('initmodel')
           call starproj%export_init3Dmodel(spproj,  params%starfile)
        case('refine3d')
           call starproj%export_class3D(spproj,  params%starfile)
        case('post')
           call starproj%export_shiny3D(spproj,  params%starfile)
        case('all')
           call starproj%export_all(spproj,  params%starfile)
        case default
           call simple_stop( " exportstar_project must have a valid export_type. Accepted values are micrograph|select|extract|class2d|initmodel|refine3d|pos or all.")
       end select


       noris = nlines(params%oritab)
       call os%new(noris)
       call os%read(params%oritab)
       if( file_exists(params%projfile) )then
           call spproj%read(params%projfile)
       endif
       call spproj%set_sp_oris(params%oritype, os)
       ! Convert Simple's single-particle project file to STAR
       call spproj%write(params%starfile)

       call spproj%kill
       call simple_end('**** EXPORTSTAR_PROJECT NORMAL STOP ****')

    end subroutine exec_exportstar_project

    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_importstar_project( self, cline )
        class(importstar_project_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(parameters) :: params
        type(sp_project) :: spproj
        character(len=LONGSTRLEN),allocatable :: starfiles(:)
        integer :: nStarfiles
        call params%new(cline)
        if( params%starfile == 'NONE')&
            call simple_stop( " importstar_project must have input file, starfile=<filename>")
        if(dir_exists(params%starfile))then
           call simple_list_files(trim(params%starfile)//'*.star', starfiles)
           nStarfiles = size(starfiles)
        else if(file_exists(params%starfile))then
           allocate(starfiles(1))
           starfiles(1)=trim(params%starfile)
           nStarfiles=1
        else
           call simple_stop( " importstar_project must have a valid input file, starfile=<filename|directory>")
        endif
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
