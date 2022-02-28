module simple_commander_starproject
include 'simple_lib.f08'
use simple_commander_base, only: commander_base
use simple_cmdline,        only: cmdline
use simple_sp_project,     only: sp_project
use simple_starproject,    only: star_project
use simple_oris,           only: oris
use simple_binoris_io,     only: binread_nlines, binread_oritab
use simple_parameters,     only: parameters, params_glob
implicit none

public :: import_starproject_commander
public :: export_starproject_commander

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

contains

  subroutine exec_import_starproject( self, cline )
  
    class(import_starproject_commander),    intent(inout)   :: self
    class(cmdline),                         intent(inout)   :: cline
    type(star_project)                                      :: starproject
    type(parameters)                                        :: params
    type(sp_project)                                        :: spproj
    integer                                                 :: it
    logical                                                 :: iteration
    character(len=3)                                        :: itchar    
    
    call cline%set('mkdir', 'yes')
    call cline%set('projname', 'project')
    
    call params%new(cline)
    
    iteration = .false.
    
    do it=999, 1, -1
        
        write(itchar, "(I0.3)") it
        
        
        if(file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
    
            iteration = .true.
            exit
            
        end if
    
    end do
    
    if(file_exists(cline%get_carg("import_dir") // "/" // "run_data.star") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_class001.mrc")) then
    
        call starproject%import_ptcls3D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "run_data.star" )
   
    else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_class001.mrc") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
        
        call starproject%import_ptcls3D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
        
    else if(iteration .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_classes.mrcs") .and. file_exists(cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star")) then
        
        call starproject%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "run_it" // itchar // "_data.star" )
    
    else if(file_exists(cline%get_carg("import_dir") // "/" // "particles.star")) then

        call starproject%import_ptcls2D(cline, spproj, cline%get_carg("import_dir") // "/" // "particles.star")

    else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")) then
        
        call starproject%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs_ctf.star")
        
    else if(file_exists(cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")) then
        
        call starproject%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "corrected_micrographs.star")
        		
    else if(file_exists(cline%get_carg("import_dir") // "/" // "micrographs.star")) then

        call starproject%import_mics(cline, spproj, cline%get_carg("import_dir") // "/" // "micrographs.star")
		
    end if
    
    call spproj%update_projinfo(cline)
    call spproj%update_compenv(cline)
    
    call spproj%write()
    
    call spproj%kill
    
    call simple_end('**** import_relion NORMAL STOP ****')
    
  end subroutine exec_import_starproject
	
  subroutine exec_export_starproject( self, cline )
  
    class(export_starproject_commander), intent(inout) :: self
    class(cmdline), intent(inout) :: cline
    type(star_project)  :: starproject
    type(parameters)     :: params
    type(sp_project)     :: spproj
    
    call params%new(cline)
    
    if( file_exists(params%projfile) )then
      call spproj%read(params%projfile)
    endif
    
   ! call starproject%create(cline, spproj)
      !  type(relion_project) :: relionproj
      !  if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'yes')
      !  if( .not. cline%defined('tiltgroups') ) call cline%set('tiltgroups', 'no')
      !  if( .not. cline%defined('reliongroups') ) call cline%set('reliongroups', 'no')
      !  if( .not. cline%defined('tiltgroupmax') ) call cline%set('tiltgroupmax', '0')
      !  if( .not. cline%defined('tiltcount') ) call cline%set('tiltcount', '0')
      !  if( .not. cline%defined('xmlloc') ) call cline%set('xmlloc', '')
      !  call params%new(cline)
      !  if( file_exists(params%projfile) )then
      !      call spproj%read(params%projfile)
      !  endif
      !  if( cline%get_rarg('reliongroups_count') .eq. 0.0) call cline%set('reliongroups_count', real(spproj%os_mic%get_noris()))
      !  call relionproj%create(cline, spproj)
   ! call starproject%kill
   
    call spproj%kill
    
    call simple_end('**** export_starproject NORMAL STOP ****')
    
  end subroutine exec_export_starproject

end module simple_commander_starproject
