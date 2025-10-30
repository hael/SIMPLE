module simple_commanders_starproject
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_starproject,    only: starproject
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters, params_glob
use simple_jiffys,         only: simple_end
use simple_nice
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_import_starproject
  contains
    procedure :: execute      => exec_import_starproject
end type commander_import_starproject

type, extends(commander_base) :: commander_export_starproject
  contains
    procedure :: execute      => exec_export_starproject
end type commander_export_starproject

type, extends(commander_base) :: commander_assign_optics_groups
  contains
    procedure :: execute      => exec_assign_optics_groups
end type commander_assign_optics_groups

contains

    subroutine exec_import_starproject( self, cline )
        class(commander_import_starproject), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(starproject)         :: starproj
        type(parameters)          :: params
        type(sp_project)          :: spproj
        integer                   :: it
        logical                   :: iteration
        character(len=3)          :: itchar
        character(len=LONGSTRLEN) :: cwd
        starproj%automode = .false.
        !show output. Defaults to false for streaming
        call starproj%set_verbose
        if(.not. cline%defined("mkdir"))       call cline%set('mkdir',       'yes' )
        if(.not. cline%defined("import_type")) call cline%set('import_type', 'auto')
        if(cline%defined("import_type") .and. index(cline%get_carg("import_type"), 'auto') > 0 ) starproj%automode = .true.
        if(starproj%automode .and. cline%defined("starfile"))             THROW_HARD('starfile should not be set with import_type=auto')
        if(.not. starproj%automode .and. .not. cline%defined("starfile")) THROW_HARD('starfile needs to be declared when not using import_type=auto')
        if(.not. dir_exists(trim(adjustl(cline%get_carg("import_dir"))))) THROW_HARD('import_dir does not exist ' // trim(adjustl(cline%get_carg("import_dir"))))
        if(cline%defined("starfile") .and. .not. file_exists(trim(adjustl(cline%get_carg("starfile"))))) THROW_HARD('starfile does not exist')
        call params%new(cline)
        if( .not. file_exists(params%projfile) ) then
            params%projfile = "workspace.simple"
            call cline%set('projfile', 'workspace.simple')
        end if
        if(file_exists(params%projfile)) then
            call spproj%read(params%projfile)
            if( trim(params%import_type).eq.'mic' )then
                if( spproj%get_nintgs() > 0 .or. spproj%get_nmovies() > 0 )then
                    THROW_HARD('The mic field of the destination PROJFILE should be empty!')
                endif
            else if( str_has_substr(trim(params%import_type),'ptcl') ) then
                 if( spproj%get_nptcls() > 0 )then
                    THROW_HARD('The ptcl fileds of the destination PROJFILE should be empty!')
                endif
            endif
        else
            THROW_HARD('Inputted projfile: '//trim(params%projfile)//' does not exist!')
        endif 
        !relative path -> absolute
        if(.not. index(cline%get_carg("import_dir"), "/") == 1) then
            call simple_getcwd(cwd)
            call cline%set('import_dir',  trim(adjustl(stemname(cwd))) // "/" // trim(adjustl(cline%get_carg("import_dir"))))
        end if
        write(logfhandle,*) ''
        if(starproj%automode) then
            write(logfhandle,*) char(9), 'using auto mode'
        else
            write(logfhandle,*) char(9), 'using manual mode'
        end if
        write(logfhandle,*) ''
        write(logfhandle,*) char(9), 'importing from ', trim(adjustl(cline%get_carg("import_dir")))
        write(logfhandle,*) ''

        if(starproj%automode) then
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
                call starproj%check_stk_params(spproj)
            else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_class001.mrc") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
                call starproj%import_ptcls3D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
                call spproj%os_ptcl2D%copy(spproj%os_ptcl3D, is_ptcl=.true.)
                call spproj%os_ptcl2D%set_all2single("e1", 0.0)
                call spproj%os_ptcl2D%set_all2single("e2", 0.0)
                call spproj%os_ptcl2D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
            else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_classes.mrcs") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
                call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
                call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
                call spproj%os_ptcl3D%set_all2single("e1", 0.0)
                call spproj%os_ptcl3D%set_all2single("e2", 0.0)
                call spproj%os_ptcl3D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
                if(file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_model.star")) then
                    call starproj%import_cls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_model.star" )
                end if
            else if(file_exists(cline%get_carg("import_dir") // "/" // "particles.star")) then
                call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "particles.star")
                call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
                call spproj%os_ptcl3D%set_all2single("e1", 0.0)
                call spproj%os_ptcl3D%set_all2single("e2", 0.0)
                call spproj%os_ptcl3D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
            else if(file_exists(cline%get_carg("import_dir") // "/" // "particles2D.star")) then
                call starproj%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "particles2D.star")
                call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
                call spproj%os_ptcl3D%set_all2single("e1", 0.0)
                call spproj%os_ptcl3D%set_all2single("e2", 0.0)
                call spproj%os_ptcl3D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
            else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")) then
                call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")
            else if(file_exists(cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")) then
                call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")
            else if(file_exists(cline%get_carg("import_dir") // "/" // "manual_pick.star")) then
                call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "manual_pick.star")
            else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs.star")) then
                call starproj%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs.star")
            else
                write(logfhandle,*) char(9), 'unable to automatically import data from folder. use manual mode'
            end if
        else
            if(index(cline%get_carg("import_type"), 'ptcl3D') > 0) then
                call starproj%import_ptcls3D(cline, spproj, trim(adjustl(cline%get_carg("starfile"))) )
                call spproj%os_ptcl2D%copy(spproj%os_ptcl3D, is_ptcl=.true.)
                call spproj%os_ptcl2D%set_all2single("e1", 0.0)
                call spproj%os_ptcl2D%set_all2single("e2", 0.0)
                call spproj%os_ptcl2D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
            else if (index(cline%get_carg("import_type"), 'ptcl2D') > 0) then
                call starproj%import_ptcls2D(cline, spproj, trim(adjustl(cline%get_carg("starfile"))) )
                call spproj%os_ptcl3D%copy(spproj%os_ptcl2D, is_ptcl=.true.)
                call spproj%os_ptcl3D%set_all2single("e1", 0.0)
                call spproj%os_ptcl3D%set_all2single("e2", 0.0)
                call spproj%os_ptcl3D%set_all2single("e3", 0.0)
                call starproj%check_stk_params(spproj)
            else if (index(cline%get_carg("import_type"), 'mic') > 0) then
                call starproj%import_mics(cline, spproj, trim(adjustl(cline%get_carg("starfile"))) )
            end if
        end if
        if(allocated(starproj%starfile%stkstates)) call spproj%report_state2stk(starproj%starfile%stkstates)
        if(spproj%get_nstks() > 0) call spproj%prune_particles()
        call spproj%update_projinfo(cline)
        call spproj%update_compenv(cline)
        call spproj%write(basename(params%projfile))
        call spproj%kill
        call starproj%kill
        call simple_end('**** IMPORT_STARPROJECT NORMAL STOP ****')
    end subroutine exec_import_starproject

    ! subroutine exec_export_starproject( self, cline )
    !     class(commander_export_starproject), intent(inout) :: self
    !     class(cmdline),                      intent(inout) :: cline
    !     type(starproject) :: starproj
    !     type(parameters)  :: params
    !     type(sp_project)  :: spproj
    !     !show output. Defaults to false for streaming
    !     call starproj%set_verbose
    !     if(.not. cline%defined("mkdir")) then
    !         call cline%set('mkdir', 'yes')
    !     end if
    !     if(.not. cline%defined("beamtilt")) then
    !         call cline%set('beamtilt', 'yes')
    !     end if
    !     call params%new(cline)
    !     call spproj%read(params%projfile)
    !     if (spproj%os_optics%get_noris() == 0) then
    !         write(logfhandle,*) ''
    !         write(logfhandle,*) char(9), "no optics groups are set in the project file. Auto assigning optics groups. You may wish to run assign_optics_groups prior to export_starproject"
    !         write(logfhandle,*) ''
    !         if(cline%get_rarg("tilt_thres") == 0) then
    !             call cline%set("tilt_thres", 0.05)
    !         end if
    !         call starproj%assign_optics(cline, spproj)
    !         call spproj%write(basename(params%projfile))
    !     end if
    !     if (spproj%os_mic%get_noris() > 0) then
    !         call starproj%export_mics(spproj)
    !     end if
    !     if (spproj%os_cls2D%get_noris() > 0) then
    !         call starproj%export_cls2D(spproj)
    !     end if
    !     if (spproj%os_ptcl2D%get_noris() > 0) then
    !         call starproj%export_ptcls2D(spproj)
    !     end if
    !     call spproj%kill
    !     call starproj%kill
    !     call simple_end('**** EXPORT_STARPROJECT NORMAL STOP ****')
    ! end subroutine exec_export_starproject

    subroutine exec_export_starproject( self, cline )
        class(commander_export_starproject), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        type(simple_nice_communicator) :: nice_communicator
        type(parameters)               :: params
        type(sp_project)               :: spproj
        if(.not. cline%defined("mkdir")) then
            call cline%set('mkdir', 'yes')
        end if
        call params%new(cline)
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        call spproj%read(params%projfile)
        call spproj%write_mics_star()
        call spproj%write_ptcl2D_star()
        call spproj%kill
        call nice_communicator%terminate()
        call simple_end('**** EXPORT_STARPROJECT NORMAL STOP ****')
    end subroutine exec_export_starproject

    subroutine exec_assign_optics_groups( self, cline )
        class(commander_assign_optics_groups), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(starproject) :: starproj
        type(parameters)  :: params
        type(sp_project)  :: spproj
        !show output. Defaults to false for streaming
        call starproj%set_verbose
        if(.not. cline%defined("mkdir")) then
            call cline%set('mkdir', 'yes')
        end if
        if(.not. cline%defined("beamtilt")) then
            call cline%set('beamtilt', 'yes')
        end if
        call params%new(cline)
        if(cline%get_rarg("tilt_thres") == 0) then
            call cline%set("tilt_thres", 0.05)
        end if
        call spproj%read(params%projfile)
        call starproj%assign_optics(cline, spproj)
        call spproj%write(basename(params%projfile))
        call spproj%kill
        call starproj%kill
        call simple_end('**** ASSIGN_OPTICS_GROUPS NORMAL STOP ****')
    end subroutine exec_assign_optics_groups

end module simple_commanders_starproject
