! set of functions to read/write/retrieve default computing parameters
module simple_computils
include 'simple_lib.f08'
#include "simple_local_flags.inc"
use json_module

implicit none

public :: get_preprocess_compconfig, get_pick_compconfig
public :: write_computing_config
private

character(len=28), parameter :: COMPDOC      = 'simple_computing_config.json'
character(len=10), parameter :: PREPROC_TMPL = 'preprocess'
character(len=4),  parameter :: PICK_TMPL    = 'pick'
character(len=17), parameter :: PREPROC_PATH = PREPROC_TMPL//'.'
character(len=11), parameter :: PICK_PATH    = PICK_TMPL//'.'
character(len=9),  parameter :: PARTITIONSTR = 'partition'
character(len=6),  parameter :: NPARTSSTR    = 'nparts'
character(len=8),  parameter :: NTHRSTR      = 'nthreads'
character(len=7),  parameter :: STREAMSTR    = 'stream_'

  contains

    ! PUBLIC ROUTINES

    ! Retrieve computing parameters for preprocessing tasks (motion correction, ctf estimation) 
    subroutine get_preprocess_compconfig( partition, nparts, nthreads, stream )
        character(len=:), allocatable, intent(inout) :: partition
        integer,                       intent(inout) :: nparts, nthreads
        logical,             optional, intent(in)    :: stream
        call get_task_compconfig(PREPROC_TMPL, partition, nparts, nthreads, stream)
    end subroutine get_preprocess_compconfig

    ! Retrieve computing parameters for picking 
    subroutine get_pick_compconfig( partition, nparts, nthreads, stream )
        character(len=:), allocatable, intent(inout) :: partition
        integer,                       intent(inout) :: nparts, nthreads
        logical,             optional, intent(in)    :: stream
        call get_task_compconfig(PICK_TMPL, partition, nparts, nthreads, stream)
    end subroutine get_pick_compconfig

    ! PRIVATE UTILITIES

    ! Retrieve computing parameters for a specific task 
    subroutine get_task_compconfig( task, partition, nparts, nthreads, stream )
        character(len=*),              intent(in)    :: task
        character(len=:), allocatable, intent(inout) :: partition
        integer,                       intent(inout) :: nparts, nthreads
        logical,             optional, intent(in)    :: stream
        type(json_file)               :: f
        character(len=:), allocatable :: path, path_prefix
        logical :: status
        ! partial path
        select case(task)
        case(PREPROC_TMPL,PICK_TMPL)
            path_prefix = trim(task)//'.'
        case DEFAULT
            THROW_HARD('Unknown task: '//task)
        end select
        if( present(stream) )then
            if( stream ) path_prefix = path_prefix//STREAMSTR
        endif
        ! read file
        call f%initialize()
        call f%load_file(filename=COMPDOC)
        if( f%failed() ) call f%check_for_errors(status)
        ! full path & retrieves
        path = path_prefix//PARTITIONSTR
        call get_str(f, path, partition)
        path = path_prefix//NPARTSSTR
        call get_int(f, path,    nparts)
        path = path_prefix//NTHRSTR
        call get_int(f, path,  nthreads)
        ! cleanup
        call f%destroy()
    end subroutine get_task_compconfig

    ! retrieves integer from json file object, returns 0 when not found
    subroutine get_int( f, path, value )
        class(json_file), intent(inout) :: f
        character(len=*), intent(in)    :: path
        integer,          intent(out)   :: value
        logical :: found
        call f%json_file_get_integer(trim(path), value, found=found)
        if( .not.found ) value = 0
    end subroutine get_int

    ! retrieves path from json file object, returns 'nil' when not found
    subroutine get_str( f, path, value )
        class(json_file),              intent(inout) :: f
        character(len=*),              intent(in)    :: path
        character(len=:), allocatable, intent(out)   :: value
        logical :: found
        call f%json_file_get_string(trim(path), value, found=found)
        if( .not.found ) value = trim(NIL)
    end subroutine get_str

    ! TESTING

    ! test routine: writes an exemplar config file
    subroutine write_computing_config()
        type(json_core)           :: json
        type(json_value), pointer :: field, entries, list
        ! init
        call json%initialize()
        call json%create_object(list,'')
        ! preprocessing
        call json%create_object(entries,PREPROC_TMPL)
        call json%add(entries, PARTITIONSTR,            'a')
        call json%add(entries, NPARTSSTR,                12)
        call json%add(entries, NTHRSTR,                   4)
        call json%add(entries, STREAMSTR//PARTITIONSTR, 'b')
        call json%add(entries, STREAMSTR//NPARTSSTR,      8)
        call json%add(entries, STREAMSTR//NTHRSTR,        6)
        call json%add(list, entries)
        ! picking 
        call json%create_object(field,  '')
        call json%create_object(entries,PICK_TMPL)
        call json%add(entries, PARTITIONSTR,            'c')
        call json%add(entries, NPARTSSTR,                 5)
        call json%add(entries, NTHRSTR,                   2)
        call json%add(entries, STREAMSTR//PARTITIONSTR, 'd')
        call json%add(entries, STREAMSTR//NPARTSSTR,      7)
        call json%add(entries, STREAMSTR//NTHRSTR,        4)
        call json%add(list, entries)
        ! write
        call json%print(list, COMPDOC)
        if( json%failed() )then
            THROW_HARD('json input/output error in write_computing_config')
        endif
        ! clean
        call json%destroy()
        nullify(list,entries)
    end subroutine write_computing_config

end module simple_computils