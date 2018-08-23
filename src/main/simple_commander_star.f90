! concrete commander: Support operations on projects for STAR-formatted files
module simple_commander_star
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_star,           only: star_project
use simple_oris,           only: oris
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters, params_glob
implicit none

public :: export_star_project_commander
public :: import_star_project_commander
public :: print_star_project_info_commander
private

type, extends(commander_base) :: export_star_project_commander
  contains
    procedure :: execute      => exec_export_star_project
end type export_star_project_commander
type, extends(commander_base) :: import_star_project_commander
  contains
    procedure :: execute      => exec_import_star_project
end type import_star_project_commander
type, extends(commander_base) ::  print_star_project_info_commander
  contains
    procedure :: execute      => exec_print_star_project_info
end type print_star_project_info_commander
#include "simple_local_flags.inc"
contains

    !> convert text (.txt) oris doc to binary (.simple)
    subroutine exec_export_star_project( self, cline )
        class(export_star_project_commander), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
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
                THROW_HARD("exportstar_project must be within a valid SIMPLE project. Projfile failed.")
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

        if( params%startype  == 'NONE')then
           write (*,*) " exportstar_project must have a valid startype."
           write (*,*) " Accepted values are micrographs|mcmicrographs|ctf_estimation|select|extract|class2d|init3dmodel|refine3d|post or all."
           write (*,*) "    movies:         for raw micrographs"
           write (*,*) "    micrographs:    for non-dose weighted micrographs"
           write (*,*) "    mcmicrographs:  for motion-corrected, dose weighted micrographs"
           write (*,*) "    ctf_estimation: for ctf estimation and micrographs"

           write (*,*) "    class2d:        for 2D class averages"
           write (*,*) "    select:         for selected 2D class averages"

           write (*,*) "    init3dmodel:    for initial particles"
           write (*,*) "    extract:        for extract dose-weighted particles"
           write (*,*) "    refine3d:       for refined 3D particles"
           write (*,*) "    post:           for post-processing"
           THROW_HARD("exportstar_project failed due to incorrect startype arg")
        end if
        !! Prepare project
!        call starproj%prepare(spproj,  params, this_starfile)

        ! Use export type to define export output variables
        select case(params%startype)
        case('movies')
           call  starproj%export_micrographs(spproj,  params%starfile)
        case('micrographs')
           call  starproj%export_micrographs(spproj,  params%starfile)
        case('mcmicrographs')
           call  starproj%export_motion_corrected_micrographs(spproj,  params%starfile)
        case('ctf_estimation')
           call starproj%export_ctf_estimation(spproj, params%starfile)
        case('select')
            call starproj%export_class2D_select(spproj,  params%starfile)
        case('extract')
           call starproj%export_extract_doseweightedptcls(spproj,  params%starfile)
        case('class2d')
           call starproj%export_class2D(spproj,  params%starfile)
        case('init3dmodel')
           call starproj%export_init3Dmodel(spproj,  params%starfile)
        case('refine3d')
           call starproj%export_class3D(spproj,  params%starfile)
        case('post')
           call starproj%export_shiny3D(spproj,  params%starfile)
        case('all')
           call starproj%export_all(spproj,  params%starfile)
        case default
           THROW_HARD("exportstar_project must have a valid startype. Accepted values are micrograph|select|extract|class2d|initmodel|refine3d|pos|all.")
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

   end subroutine exec_export_star_project

    !> for importing STAR formatted project files to *.simple
   subroutine exec_import_star_project( self, cline )
       use simple_oris,      only: oris
       use simple_nrtxtfile, only: nrtxtfile
       use simple_binoris_io ! use all in there
       class(import_star_project_commander), intent(inout) :: self
       class(cmdline),                       intent(inout) :: cline
       type(parameters) :: params
       type(sp_project) :: spproj
       type(star_project):: starproj
       type(ctfparams)  :: ctfvars
       type(oris)       :: os
       character(len=LONGSTRLEN),allocatable :: line(:)
       character(len=LONGSTRLEN),allocatable :: starfiles(:)
       character(len=:),      allocatable :: phaseplate, boxf_abspath
       character(len=LONGSTRLEN), allocatable :: boxfnames(:)
       integer :: nStarfiles,istar,i, ndatlines,nrecs, nboxf, nmovf
       logical :: l_starfile, inputted_boxtab ,inputted_smpd, inputted_startype
       call params%new(cline)

       ! parameter input management
       l_starfile = cline%defined('starfile')
       if( .not. l_starfile)then
           THROW_HARD('Star project file argument empty or not set, e.g. starfile=<filename> or <directory>')
       end if
       if(dir_exists(params%starfile))then
           call simple_list_files(trim(params%starfile)//'*.star', starfiles)
           nStarfiles = size(starfiles)
           write(*,*) " Importing star project :", nStarfiles, " found "
       else if(file_exists(params%starfile))then
           allocate(starfiles(1))
           starfiles(1)=trim(params%starfile)
           nStarfiles=1
       else
           write(*,*) " "
           THROW_HARD('Importing star project must have a valid input file, starfile=<filename|directory>')
       endif
       inputted_startype = cline%defined('startype')
       if(inputted_startype)then
           write (*,*) "Testing for valid startype:", trim(params%startype)
           if( RE_match(params%startype,&
               '(m|movies|micrographs|'&
               'ctf|ctf_estimation|ctfparams|'//&
               'p|ptcl|particles|'//&
               'cavg|classaverages)') /=0 )then
               inputted_startype=.false.
               write (*,*) "ERROR: Invalid startype:", trim(params%startype)
           endif
       endif
       if( .not. inputted_startype)then
           write (*,*) " import_starproject valid startypes:"
           write (*,*) " Accepted values are m|movies|micrographs ==> import micrographs"
           write (*,*) "                     ctf|ctf_estimation   ==> import mic + ctf params"
           write (*,*) "                     p|ptcl|particles     ==> import particles "
           write (*,*) "                     cavgs|classaverages  ==> import class averages "

           THROW_HARD('import_starproject; `startype` argument empty or not set.')
       end if
       inputted_smpd = cline%defined('smpd')
       if(.not. inputted_smpd)then
           THROW_HARD('Importing star project must have `smpd` as argument.')
       endif

       inputted_boxtab = cline%defined('boxtab')

       ! project file management
       if( .not. file_exists(trim(params%projfile)) )then
           write(*,*) 'Project file: ', trim(params%projfile), ' does not exist!'
           write(*,*)
           THROW_HARD("exec_import_starproject; not in SIMPLE project dir. ")
       endif
       !! Read existing SIMPLE project
       call spproj%read(params%projfile)

       ! Prepare STAR project module
       call starproj%prepareimport( spproj, params, starfiles(1))
       ndatlines = starproj%get_ndatalines()
       nrecs     = starproj%get_nrecs_per_line()
       call starproj%check_temp_files("importstar commander")

       !! Check number of records on data line to validate starfile
       if(nrecs == 0)then
           THROW_HARD("exec_import_starproject; Unable to read header records in STAR file.")
       end if
       if(ndatlines == 0)then
           THROW_HARD("exec_import_starproject; Unable to read data line records in STAR file.")
       end if


       !! vvvvv TODO vvvvv

       select case(params%startype)
       case('m')
           call  starproj%import_micrographs(spproj, params, cline, params%starfile)
       case('movies')
           call  starproj%import_micrographs(spproj, params, cline, params%starfile)
       case('micrographs')
           call  starproj%import_micrographs(spproj, params, cline, params%starfile)
  !     case('mcmicrographs')
  !         call  starproj%import_motion_corrected_micrographs(spproj, params, cline, params%starfile)
       case('ctf')
           call starproj%import_ctf_estimation(spproj, params, cline, params%starfile)
       case('ctf_estimation')
           call starproj%import_ctf_estimation(spproj, params, cline, params%starfile)
  !     case('select')
  !!         call starproj%import_class2D_select(spproj, params, cline,params%starfile)
  !     case('extract')
  !         call starproj%import_extract_doseweightedptcls(spproj,params, cline, params%starfile)
       case('cavgs')
          call starproj%import_cavgs(spproj, params, cline, params%starfile)
       case('classaverages')
          call starproj%import_cavgs(spproj, params, cline, params%starfile)
  !     case('init3dmodel')
  !         call starproj%import_init3Dmodel(spproj, params, cline, params%starfile)
       case('p')
           call starproj%import_particles(spproj, params, cline, params%starfile)
       case('ptcls')
           call starproj%import_particles(spproj, params, cline, params%starfile)
     case('particles')
           call starproj%import_particles(spproj, params, cline, params%starfile)
  !     case('post')
  !         call starproj%import_shiny3D(spproj, params, cline, params%starfile)
  !     case('all')
  !         call starproj%import_all(spproj, params, cline, params%starfile)
       case default
           THROW_HARD("import_starproject must have a valid startype. ")
       end select

       call starproj%print_info

       ! Import STAR filename
       ! do istar=1, nStarfiles
       !      call starproj%read(starfiles(i), spproj, params)
       ! end do

       ! Copy to SIMPLE's project format

       call spproj%write
       call spproj%kill
       call simple_end('**** import_starproject NORMAL STOP ****')
   end subroutine exec_import_star_project


    !> convert binary (.simple) oris doc to text (.txt)
    subroutine exec_print_star_project_info( self, cline )
        class(print_star_project_info_commander), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
         type(parameters) :: params
         type(star_project) :: starproj
         call params%new(cline)
         print *," Reading star-formatted file: ", trim(params%starfile)
         call starproj%read(params%starfile)
         call starproj%print_info
         call starproj%kill
         call simple_end('**** PRINT_PROJECT_STAR_INFO NORMAL STOP ****')
    end subroutine exec_print_star_project_info



end module simple_commander_star
