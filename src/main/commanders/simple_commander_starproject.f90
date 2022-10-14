module simple_commander_starproject
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_starproject,    only: starproject
use simple_oris,           only: oris
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters, params_glob
use simple_syslib,         only: simple_getcwd
implicit none

public :: import_starproject_commander
public :: export_starproject_commander
public :: assign_optics_groups_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: import_starproject_commander
  contains
    procedure :: execute      => exec_import_starproject
end type import_starproject_commander

type, extends(commander_base) :: export_starproject_commander
  contains
    procedure :: execute      => exec_export_starproject
end type export_starproject_commander

type, extends(commander_base) :: assign_optics_groups_commander
  contains
    procedure :: execute      => exec_assign_optics_groups
end type assign_optics_groups_commander

contains

    subroutine exec_import_starproject( self, cline )
        class(import_starproject_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(starproject)         :: starproj
        type(parameters)          :: params
        type(sp_project)          :: spproj
        integer                   :: it
        logical                   :: iteration
        character(len=3)          :: itchar    
        character(len=LONGSTRLEN) :: cwd
        !show output. Defaults to false for streaming
        call starproj%set_verbose
        if(.not. cline%defined("mkdir")) then
            call cline%set('mkdir', 'yes')
        end if
        call params%new(cline)
        call spproj%read(params%projfile)
        if(.not. index(cline%get_carg("import_dir"), "/") == 1) then
            call simple_getcwd(cwd)
            call cline%set('import_dir',  trim(adjustl(stemname(cwd))) // "/" // trim(adjustl(cline%get_carg("import_dir"))))
        end if
        if(dir_exists(trim(adjustl(cline%get_carg("import_dir"))))) then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), 'importing from ', trim(adjustl(cline%get_carg("import_dir")))
            write(logfhandle,*) ''
        else
            THROW_HARD('folder does not exist ' // trim(adjustl(cline%get_carg("import_dir"))))
        end if
        iteration = .false.
        do it=999, 1, -1
            write(itchar, "(I0.3)") it
            if(file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
                iteration = .true.
                exit
            end if
        end do
        if(file_exists(cline%get_carg("import_dir") // "/" // "run_data.star") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_class001.mrc")) then
            call starproj%import_ptcls3D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "run_data.star" )
            call spproj%os_ptcl2D%copy(spproj%os_ptcl3D, is_ptcl=.true.)
            call spproj%os_ptcl2D%set_all2single("e1", 0.0)
            call spproj%os_ptcl2D%set_all2single("e2", 0.0)
            call spproj%os_ptcl2D%set_all2single("e3", 0.0)
        else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_class001.mrc") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
            call starproj%import_ptcls3D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
            call spproj%os_ptcl2D%copy(spproj%os_ptcl3D, is_ptcl=.true.)
            call spproj%os_ptcl2D%set_all2single("e1", 0.0)
            call spproj%os_ptcl2D%set_all2single("e2", 0.0)
            call spproj%os_ptcl2D%set_all2single("e3", 0.0)
        else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_classes.mrcs") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
            call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
            call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
            call spproj%os_ptcl3D%set_all2single("e1", 0.0)
            call spproj%os_ptcl3D%set_all2single("e2", 0.0)
            call spproj%os_ptcl3D%set_all2single("e3", 0.0)
            if(file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_model.star")) then
                call starproj%import_cls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_model.star" )
            end if
        else if(file_exists(cline%get_carg("import_dir") // "/" // "particles.star")) then
            call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "particles.star")
            call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
            call spproj%os_ptcl3D%set_all2single("e1", 0.0)
            call spproj%os_ptcl3D%set_all2single("e2", 0.0)
            call spproj%os_ptcl3D%set_all2single("e3", 0.0)
        else if(file_exists(cline%get_carg("import_dir") // "/" // "particles2D.star")) then
            call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "particles2D.star")
            call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
            call spproj%os_ptcl3D%set_all2single("e1", 0.0)
            call spproj%os_ptcl3D%set_all2single("e2", 0.0)
            call spproj%os_ptcl3D%set_all2single("e3", 0.0)
        else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")) then
            call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")
        else if(file_exists(cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")) then
            call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")
        else if(file_exists(cline%get_carg("import_dir") // "/" // "manual_pick.star")) then
            call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "manual_pick.star")
        else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs.star")) then
            call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs.star")
        end if
        call spproj%update_projinfo(cline)
        call spproj%update_compenv(cline)
        call spproj%write(basename(params%projfile))
        call spproj%kill
        call simple_end('**** IMPORT_STARPROJECT NORMAL STOP ****')
    end subroutine exec_import_starproject
	
    subroutine exec_export_starproject( self, cline )
        class(export_starproject_commander), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(starproject) :: starproj
        type(parameters)  :: params
        type(sp_project)  :: spproj
        !show output. Defaults to false for streaming
        call starproj%set_verbose
        if(.not. cline%defined("mkdir")) then
            call cline%set('mkdir', 'yes')
        end if
        call params%new(cline)
        call spproj%read(params%projfile)
        if (spproj%os_optics%get_noris() == 0) then
            write(logfhandle,*) ''
            write(logfhandle,*) char(9), "no optics groups are set in the project file. Auto assigning optics groups. You may wish to run assign_optics_groups prior to export_starproject"
            write(logfhandle,*) ''
            if(cline%get_rarg("tilt_thres") == 0) then
                call cline%set("tilt_thres", 0.05)
            end if
            call starproj%assign_optics(cline, spproj)
            call spproj%write(basename(params%projfile))
        end if
        if (spproj%os_mic%get_noris() > 0) then
            call starproj%export_mics(cline, spproj)
        end if
        if (spproj%os_cls2D%get_noris() > 0) then
            call starproj%export_cls2D(spproj)
        end if
        if (spproj%os_ptcl2D%get_noris() > 0) then
            call starproj%export_ptcls2D(cline, spproj)
        end if
        call spproj%kill
        call simple_end('**** EXPORT_STARPROJECT NORMAL STOP ****')
    end subroutine exec_export_starproject

    subroutine exec_assign_optics_groups( self, cline )
        class(assign_optics_groups_commander), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(starproject) :: starproj
        type(parameters)  :: params
        type(sp_project)  :: spproj
        !show output. Defaults to false for streaming
        call starproj%set_verbose
        call cline%set('mkdir', 'yes')
        call params%new(cline)
        if(cline%get_rarg("tilt_thres") == 0) then
            call cline%set("tilt_thres", 0.05)
        end if
        call spproj%read(params%projfile)
        call starproj%assign_optics(cline, spproj)
        call spproj%write(basename(params%projfile))
        call spproj%kill
        call simple_end('**** ASSIGN_OPTICS_GROUPS NORMAL STOP ****')
    end subroutine exec_assign_optics_groups

end module simple_commander_starproject
