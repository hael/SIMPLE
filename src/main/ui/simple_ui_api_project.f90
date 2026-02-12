!@descr: "project" UI api (concrete implementation)
module simple_ui_api_project
use simple_ui_program, only: ui_program
use simple_ui_hash,    only: ui_hash
use simple_ui_utils,   only: add_ui_program
implicit none

type(ui_program), target :: export_relion
type(ui_program), target :: export_starproject
type(ui_program), target :: extract_subproj
type(ui_program), target :: import_boxes
type(ui_program), target :: import_cavgs
type(ui_program), target :: import_movies
type(ui_program), target :: import_particles
type(ui_program), target :: import_starproject
type(ui_program), target :: merge_projects
type(ui_program), target :: new_project
type(ui_program), target :: print_project_field
type(ui_program), target :: print_project_info
type(ui_program), target :: prune_project
type(ui_program), target :: replace_project_field
type(ui_program), target :: selection
type(ui_program), target :: update_project
type(ui_program), target :: zero_project_shifts

contains

    subroutine register_ui_project(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('export_relion',         export_relion,         prgtab)
        call add_ui_program('export_starproject',    export_starproject,    prgtab)
        call add_ui_program('extract_subproj',       extract_subproj,       prgtab)
        call add_ui_program('import_boxes',          import_boxes,          prgtab)
        call add_ui_program('import_cavgs',          import_cavgs,          prgtab)
        call add_ui_program('import_movies',         import_movies,         prgtab)
        call add_ui_program('import_particles',      import_particles,      prgtab)
        call add_ui_program('import_starproject',    import_starproject,    prgtab)
        call add_ui_program('merge_projects',        merge_projects,        prgtab)
        call add_ui_program('new_project',           new_project,           prgtab)
        call add_ui_program('print_project_field',   print_project_field,   prgtab)
        call add_ui_program('print_project_info',    print_project_info,    prgtab)
        call add_ui_program('prune_project',         prune_project,         prgtab)
        call add_ui_program('replace_project_field', replace_project_field, prgtab)
        call add_ui_program('selection',             selection,             prgtab)
        call add_ui_program('update_project',        update_project,        prgtab)
        call add_ui_program('zero_project_shifts',   zero_project_shifts,   prgtab)
    end subroutine register_ui_project

end module simple_ui_api_project
