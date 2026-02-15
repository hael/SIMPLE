!@descr: execution of project management commanders
module simple_exec_project
use simple_cmdline,                 only: cmdline
use simple_commanders_relion,       only: commander_export_relion
use simple_commanders_starproject,  only: commander_import_starproject, commander_export_starproject
use simple_commanders_project_core, only: commander_new_project, commander_update_project, commander_print_project_info,&
commander_print_project_field, commander_replace_project_field, commander_selection, commander_merge_projects,&
commander_extract_subproj
use simple_commanders_project_mov,  only: commander_import_movies, commander_write_mic_filetab
use simple_commanders_project_ptcl, only: commander_zero_project_shifts, commander_import_boxes,&
commander_import_particles, commander_prune_project_distr
use simple_commanders_project_cls,  only: commander_import_cavgs
implicit none

public :: exec_project_commander
private

type(commander_export_relion)         :: xexport_relion
type(commander_export_starproject)    :: xexport_starproject
type(commander_extract_subproj)       :: xextract_subproj
type(commander_import_boxes)          :: ximport_boxes
type(commander_import_cavgs)          :: ximport_cavgs
type(commander_import_movies)         :: ximport_movies
type(commander_import_particles)      :: ximport_particles
type(commander_import_starproject)    :: ximport_starproject
type(commander_merge_projects)        :: xmerge_projects
type(commander_new_project)           :: xnew_project
type(commander_print_project_field)   :: xprint_project_field
type(commander_print_project_info)    :: xprint_project_info
type(commander_prune_project_distr)   :: xprune_project
type(commander_replace_project_field) :: xreplace_project_field
type(commander_selection)             :: xselection
type(commander_update_project)        :: xupdate_project
type(commander_zero_project_shifts)   :: xzero_project_shifts
type(commander_write_mic_filetab)     :: xwrite_mic_filetab

contains

    subroutine exec_project_commander(which, cline, l_silent, l_did_execute)
        character(len=*),    intent(in)    :: which
        class(cmdline),      intent(inout) :: cline
        logical,             intent(out)   :: l_silent
        logical,             intent(out)   :: l_did_execute
        l_silent      = .false.
        l_did_execute = .true.
        select case(trim(adjustl(which)))
            case( 'export_relion' )
                call xexport_relion%execute(cline)
            case( 'export_starproject' )
                call xexport_starproject%execute(cline)
            case( 'extract_subproj' )
                call xextract_subproj%execute(cline)
            case( 'import_boxes' )
                call ximport_boxes%execute(cline)
            case( 'import_cavgs' )
                call ximport_cavgs%execute(cline)
            case( 'import_movies' )
                call ximport_movies%execute(cline)
            case( 'import_particles' )
                call ximport_particles%execute(cline)
            case( 'import_starproject' )
                call ximport_starproject%execute(cline)
            case( 'merge_projects' )
                call xmerge_projects%execute(cline)
            case( 'new_project' )
                call xnew_project%execute(cline)
            case( 'print_project_field' )
                call xprint_project_field%execute(cline)
                l_silent = .true.
            case( 'print_project_info' )
                call xprint_project_info%execute(cline)
                l_silent = .true.
            case( 'prune_project' )
                call xprune_project%execute(cline)
            case( 'replace_project_field' )
                call xreplace_project_field%execute(cline)
            case( 'selection', 'report_selection' )
                call xselection%execute(cline)
            case( 'update_project' )
                call xupdate_project%execute(cline)
            case( 'zero_project_shifts' )
                call xzero_project_shifts%execute(cline)
            case( 'write_mic_filetab' )
                call xwrite_mic_filetab%execute(cline)
            case default
                l_did_execute = .false.
        end select
    end subroutine exec_project_commander

end module simple_exec_project
