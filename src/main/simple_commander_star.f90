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

           call simple_stop( " exportstar_project failed due to incorrect startype arg")
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
           call starproj%export_ctf_estimation(spproj,  params%starfile)
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
           call simple_stop( " exportstar_project must have a valid startype. Accepted values are micrograph|select|extract|class2d|initmodel|refine3d|pos or all.")
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
        logical :: l_starfile, inputted_boxtab
        call params%new(cline)
        ! parameter input management
        l_starfile = cline%defined('starfile')
        if( .not. l_starfile)then
            write(*,*) 'Star project file argument empty or not set, e.g. starfile=<filename> or <directory>'
            call simple_stop( "ERROR! simple_commander_star :: exec_importstar_project ")
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
           write(*,*) " Importing star project must have a valid input file, starfile=<filename|directory>"
           call simple_stop('ERROR! simple_commander_star :: exec_importstar_project')
       endif
       inputted_boxtab = cline%defined('boxtab')

       ! project file management
       if( .not. file_exists(trim(params%projfile)) )then
           write(*,*) 'Project file: ', trim(params%projfile), ' does not exist!'
           HALT_NOW("exec_importstar_project; not in SIMPLE project dir")
       endif
        !! Read existing SIMPLE project
        call spproj%read(params%projfile)

       ! Prepare STAR project module
       call starproj%prepareimport( spproj, params, starfiles(1))
       ndatlines = starproj%get_ndatalines()
       nrecs     = starproj%get_nrecs_per_line()

       if(nrecs == 0)then
           HALT_NOW("exec_importstar_project; Unable to read header records in STAR file.")
       end if
        if(ndatlines == 0)then
           HALT_NOW("exec_importstar_project; Unable to read data line records in STAR file.")
       end if

        if( params%startype  == 'NONE')then
           write (*,*) " importstar_project valid startypes:"
           write (*,*) " Accepted values are movies|micrographs|mcmicrographs|ctf_estimation|select|extract|class2d|init3dmodel|refine3d|post or all."
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

           
        end if


!! vvvvv TODO vvvvv

        select case(params%startype)
        case('movies')
           call  starproj%import_micrographs(spproj, params, params%starfile)
        case('micrographs')
           call  starproj%import_micrographs(spproj, params, params%starfile)
        case('mcmicrographs')
           call  starproj%import_motion_corrected_micrographs(spproj, params, params%starfile)
        case('ctf_estimation')
           call starproj%import_ctf_estimation(spproj, params, params%starfile)
        case('select')
            call starproj%import_class2D_select(spproj,params,  params%starfile)
        case('extract')
           call starproj%import_extract_doseweightedptcls(spproj,params,  params%starfile)
        case('class2d')
           call starproj%import_class2D(spproj, params, params%starfile)
        case('init3dmodel')
           call starproj%import_init3Dmodel(spproj, params, params%starfile)
        case('refine3d')
           call starproj%import_class3D(spproj, params, params%starfile)
        case('post')
           call starproj%import_shiny3D(spproj, params, params%starfile)
        case('all')
           call starproj%import_all(spproj, params, params%starfile)
        case default
           call simple_stop( " importstar_project must have a valid startype. Accepted values are micrograph|select|extract|class2d|initmodel|refine3d|pos or all.")
       end select




!!******* Import movies + ctf_estimation
       if( cline%defined('phaseplate') )then
            phaseplate = cline%get_carg('phaseplate')
        else
            allocate(phaseplate, source='no')
        endif
        ctfvars%smpd  = params%smpd
        ctfvars%kv    = params%kv
        ctfvars%cs    = params%cs
        ctfvars%fraca = params%fraca
        select case(params%ctf)
            case('yes')
                ctfvars%ctfflag = 1
            case('no')
                ctfvars%ctfflag = 0
            case('flip')
                ctfvars%ctfflag = 2
            case DEFAULT
                write(*,*) 'ctf flag params%ctf: ', params%ctf
                stop 'ERROR! ctf flag not supported; commander_project :: import_movies'
        end select
        ctfvars%l_phaseplate = .false.
        if( trim(params%phaseplate) .eq. 'yes' ) ctfvars%l_phaseplate = .true.
        ! update project info
        call spproj%update_projinfo( cline )
        ! updates segment
        call spproj%add_movies(params%filetab, ctfvars)
        ! add boxtab
        if( file_exists(params%boxtab) )then
            call read_filetable(params%boxtab, boxfnames)
            nboxf = size(boxfnames)
            nmovf = nlines(params%filetab)
            if( nboxf /= nmovf )then
                write(*,*) '# boxfiles: ', nboxf
                write(*,*) '# movies  : ', nmovf
                stop 'ERROR! # boxfiles .ne. # movies; commander_project :: exec_import_movies'
            endif
            do i=1,nmovf
                call simple_abspath(boxfnames(i), boxf_abspath, errmsg='commander_project :: exec_import_movies')
                call spproj%os_mic%set(i, 'boxfile', boxf_abspath)
            end do
        end if



!! ***** Class averages
!        call spproj%add_cavgs2os_out(params%stk, params%smpd)

!! ******  Import boxes
        ! ! get boxfiles into os_mic
        ! call read_filetable(params%boxtab, boxfnames)
        ! nboxf   = size(boxfnames)
        ! nos_mic = spproj%os_mic%get_noris()
        ! if( nboxf /= nos_mic )then
        !     write(*,*) '# boxfiles       : ', nboxf
        !     write(*,*) '# os_mic entries : ', nos_mic
        !     stop 'ERROR! # boxfiles .ne. # os_mic entries; commander_project :: exec_import_boxes'
        ! endif
        ! do i=1,nos_mic
        !     call abspath(trim(boxfnames(i)), boxf_abspath, 'commander_project :: exec_import_movies')
        !     call spproj%os_mic%set(i, 'boxfile', boxf_abspath)
        ! end do

!! ******  Import  particles

        if( nrecs < 1 .or. nrecs > 4 .or. nrecs == 2 )then
           write(*,*) 'unsupported nr of rec:s in plaintexttab'
           stop 'commander_starproject :: exec_extract_ptcls'
        endif
!   call os%new(ndatlines)
! allocate( line(nrecs) )
!             do i=1,ndatlines
!                 call starproj%readNextDataLine(line)
!                 select case(params%dfunit)
!                     case( 'A' )
!                         line(1) = line(1)/1.0e4
!                         if( nrecs > 1 )  line(2) = line(2)/1.0e4
!                     case( 'microns' )
!                         ! nothing to do
!                     case DEFAULT
!                         stop 'unsupported dfunit; commander_project :: exec_extract_ptcls'
!                 end select
!                 select case(params%angastunit)
!                     case( 'radians' )
!                         if( nrecs == 3 ) line(3) = rad2deg(line(3))
!                     case( 'degrees' )
!                         ! nothing to do
!                     case DEFAULT
!                         stop 'unsupported angastunit; commander_project :: exec_extract_ptcls'
!                 end select
!                 select case(params%phshiftunit)
!                     case( 'radians' )
!                         ! nothing to do
!                     case( 'degrees' )
!                         if( nrecs == 4 ) line(4) = deg2rad(line(4))
!                     case DEFAULT
!                         stop 'unsupported phshiftunit; commander_project :: exec_extract_ptcls'
!                 end select
!  call os%set(i, 'dfx', line(1))
!                 if( nrecs > 1 )then
!                     call os%set(i, 'dfy', line(2))
!                     call os%set(i, 'angast', line(3))
!                 endif
!                 if( nrecs > 3 )then
!                     call os%set(i, 'phshift', line(4))
!                 endif
!             end do


        ! Import STAR filename
        ! do istar=1, nStarfiles
        !      call starproj%read(starfiles(i), spproj, params)
        ! end do

        ! Copy to SIMPLE's project format

        call spproj%write
        call spproj%kill
        call simple_end('**** IMPORTSTAR_PROJECT NORMAL STOP ****')
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
