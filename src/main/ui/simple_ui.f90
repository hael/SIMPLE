!@descr: the main user interface module 
module simple_ui
! core helpers
use simple_core_module_api
use simple_ansi_ctrls
use simple_ui_params_common
use simple_ui_hash,         only: ui_hash
use simple_ui_program,      only: UI_DISPLAY_NAME_MAX_LEN, UI_SUMMARY_MIN_LEN, UI_SUMMARY_MAX_LEN, &
    &category_descriptor, ui_program, ui_program_input
use simple_ui_param,        only: UI_PLACEHOLDER_MAX_LEN, ui_placeholder_is_standard
! program table grouping helpers
use simple_ui_simple_group, only: add_simple_programs
use simple_ui_stream_group, only: add_stream_programs
use simple_ui_single_group, only: add_single_programs
use simple_ui_test_group,   only: add_test_programs
implicit none

public :: make_ui, make_test_ui, get_prg_ptr, get_test_prg_ptr
public :: list_simple_prgs_in_ui, list_simple_test_prgs_in_ui, list_stream_prgs_in_ui, list_single_prgs_in_ui
public :: print_ui_json, write_ui_json
public :: print_stream_ui_json, validate_ui_json, validate_ui_presentation
private
#include "simple_local_flags.inc"

character(len=26), parameter :: UI_FNAME = 'simple_ui.json'
character(len=*),  parameter :: UI_JSON_REAL_FORMAT = '(ss,G0.6)'
logical,           parameter :: DEBUG    = .false.
type(ui_hash)                :: prgtab, tsttab
type(string), allocatable    :: prgnames(:)
type(string), allocatable    :: tstnames(:)

contains
    
    ! public class methods

    subroutine make_ui
        call set_ui_params
        ! SIMPLE PROGRAMS
        call add_simple_programs(prgtab)
        ! SIMPLE STREAM PROGRAMS
        call add_stream_programs(prgtab)
        ! SINGLE PROGRAMS
        call add_single_programs(prgtab)
        prgnames = prgtab%keys_sorted()
        call validate_category_metadata(prgtab, prgnames)
        call validate_ui_presentation
    end subroutine make_ui

    subroutine make_test_ui
        call set_ui_params
        ! SIMPLE TEST PROGRAMS
        call add_test_programs(tsttab)
        tstnames = tsttab%keys_sorted()
        call validate_category_metadata(tsttab, tstnames)
    end subroutine make_test_ui

    subroutine get_prg_ptr( which_program, ptr2prg )
        class(string),             intent(in)    :: which_program
        type(ui_program), pointer, intent(inout) :: ptr2prg
        ptr2prg => null()
        call prgtab%get_ui_program(which_program, ptr2prg)
    end subroutine get_prg_ptr

    subroutine validate_ui_presentation
        use simple_linked_list, only: linked_list, list_iterator
        type(ui_program), pointer :: ptr2prg
        integer                   :: iprg
        if( .not. allocated(prgnames) ) return
        do iprg = 1, size(prgnames)
            call prgtab%get_ui_program(prgnames(iprg), ptr2prg)
            if( len_trim(ptr2prg%summary%to_char()) < UI_SUMMARY_MIN_LEN .or. &
                &len_trim(ptr2prg%summary%to_char()) > UI_SUMMARY_MAX_LEN )then
                THROW_HARD('ui_program summary must contain 30-100 characters: '//ptr2prg%name%to_char())
            endif
            if( len_trim(ptr2prg%display_name%to_char()) == 0 .or. &
                &len_trim(ptr2prg%display_name%to_char()) > UI_DISPLAY_NAME_MAX_LEN )then
                THROW_HARD('ui_program display_name must contain 1-100 characters: '//ptr2prg%name%to_char())
            endif
            call validate_input_list(ptr2prg%img_ios)
            call validate_input_list(ptr2prg%file_ios)
            call validate_input_list(ptr2prg%parm_ios)
            call validate_input_list(ptr2prg%srch_ctrls)
            call validate_input_list(ptr2prg%filt_ctrls)
            call validate_input_list(ptr2prg%mask_ctrls)
            call validate_input_list(ptr2prg%comp_ctrls)
        enddo
    contains
        subroutine validate_input_list( list )
            type(linked_list), intent(in) :: list
            type(list_iterator)           :: iterator
            class(*), allocatable         :: item
            iterator = list%begin()
            do while( iterator%has_value() )
                call iterator%getter(item)
                select type( input => item )
                    type is( ui_program_input )
                        if( len_trim(input%param%placeholder%to_char()) > UI_PLACEHOLDER_MAX_LEN .or. &
                            &.not. ui_placeholder_is_standard(input%param%key%to_char(), input%param%keytype%to_char(), &
                            &input%param%placeholder%to_char(), input%param%units%to_char(), &
                            &input%param%example_value%to_char()) )then
                            THROW_HARD('ui_param placeholder is not standardized: '//input%param%key%to_char())
                        endif
                end select
                if( allocated(item) ) deallocate(item)
                call iterator%next()
            enddo
        end subroutine validate_input_list
    end subroutine validate_ui_presentation

    subroutine validate_category_metadata( program_table, program_names )
        class(ui_hash), intent(in) :: program_table
        type(string),   intent(in) :: program_names(:)
        type(ui_program), pointer :: first_program, second_program
        type(string) :: errmsg
        integer :: ifirst, isecond
        do ifirst = 1, size(program_names)
            call program_table%get_ui_program(program_names(ifirst), first_program)
            if( len_trim(first_program%category%to_char()) == 0 .or. &
                &len_trim(first_program%category_display_name%to_char()) == 0 .or. &
                &first_program%category_order <= 0 )then
                THROW_HARD('incomplete category metadata for program: '//first_program%name%to_char())
            endif
            do isecond = ifirst + 1, size(program_names)
                call program_table%get_ui_program(program_names(isecond), second_program)
                if( trim(first_program%executable%to_char()) /= trim(second_program%executable%to_char()) ) cycle
                if( trim(first_program%category%to_char()) /= trim(second_program%category%to_char()) )then
                    ! Distinct categories in one executable must have distinct display order.
                    if( first_program%category_order == second_program%category_order )then
                        errmsg = 'duplicate category order in executable '//trim(first_program%executable%to_char())//': '// &
                            &trim(first_program%category%to_char())//' and '//trim(second_program%category%to_char())
                        THROW_HARD(errmsg%to_char())
                    endif
                    cycle
                endif
                if( trim(first_program%category_display_name%to_char()) /= &
                    &trim(second_program%category_display_name%to_char()) .or. &
                    &first_program%category_order /= second_program%category_order )then
                    THROW_HARD('inconsistent metadata for category: '//first_program%category%to_char())
                endif
            enddo
        enddo
    end subroutine validate_category_metadata

    subroutine get_test_prg_ptr( which_program, ptr2prg )
        class(string),             intent(in)    :: which_program
        type(ui_program), pointer, intent(inout) :: ptr2prg
        ptr2prg => null()
        call tsttab%get_ui_program(which_program, ptr2prg)
    end subroutine get_test_prg_ptr

    subroutine list_simple_prgs_in_ui
        call print_registered_programs(prgtab, prgnames, 'simple_exec')
    end subroutine list_simple_prgs_in_ui

    subroutine list_simple_test_prgs_in_ui
        call print_registered_programs(tsttab, tstnames)
    end subroutine list_simple_test_prgs_in_ui

    subroutine list_stream_prgs_in_ui
        call print_registered_programs(prgtab, prgnames, 'simple_stream')
    end subroutine list_stream_prgs_in_ui

    subroutine list_single_prgs_in_ui
        call print_registered_programs(prgtab, prgnames, 'single_exec')
    end subroutine list_single_prgs_in_ui

    subroutine print_registered_programs( program_table, program_names, executable )
        class(ui_hash), intent(in) :: program_table
        type(string),   intent(in) :: program_names(:)
        character(len=*), optional, intent(in) :: executable
        type(ui_program), pointer :: program
        type(category_descriptor), allocatable :: categories(:)
        type(category_descriptor) :: category
        integer :: icategory, iprogram, jcategory
        logical :: found

        do iprogram = 1, size(program_names)
            call program_table%get_ui_program(program_names(iprogram), program)
            if( present(executable) )then
                if( trim(program%executable%to_char()) /= trim(executable) ) cycle
            endif
            category%id           = program%category%to_char()
            category%display_name = program%category_display_name%to_char()
            category%order        = program%category_order
            found = .false.
            if( allocated(categories) )then
                do icategory = 1, size(categories)
                    if( trim(categories(icategory)%id) == trim(category%id) )then
                        if( trim(categories(icategory)%display_name) /= trim(category%display_name) .or. &
                            &categories(icategory)%order /= category%order )then
                            THROW_HARD('inconsistent metadata for category: '//trim(category%id))
                        endif
                        found = .true.
                        exit
                    endif
                enddo
            endif
            if( .not. found )then
                if( allocated(categories) )then
                    categories = [categories, category]
                else
                    allocate(categories(1))
                    categories(1) = category
                endif
            endif
        enddo

        if( .not. allocated(categories) ) return
        do icategory = 1, size(categories) - 1
            do jcategory = icategory + 1, size(categories)
                if( categories(jcategory)%order < categories(icategory)%order )then
                    category = categories(icategory)
                    categories(icategory) = categories(jcategory)
                    categories(jcategory) = category
                endif
            enddo
        enddo

        do icategory = 1, size(categories)
            write(logfhandle,'(A)') format_str(trim(categories(icategory)%display_name)//':', C_UNDERLINED)
            do iprogram = 1, size(program_names)
                call program_table%get_ui_program(program_names(iprogram), program)
                if( present(executable) )then
                    if( trim(program%executable%to_char()) /= trim(executable) ) cycle
                endif
                if( trim(program%category%to_char()) == trim(categories(icategory)%id) )then
                    write(logfhandle,'(A)') program%name%to_char()
                endif
            enddo
            write(logfhandle,'(A)') ''
        enddo
    end subroutine print_registered_programs

    subroutine print_ui_json
        use json_module
        type(json_core)           :: json
        type(json_value), pointer :: all_programs, program_entry
        type(ui_program), pointer :: ptr2prg
        integer                   :: iprg
        ! JSON init
        call json%initialize(real_format=UI_JSON_REAL_FORMAT)
        ! create object of program entries
        call json%create_object(all_programs, 'SIMPLE_UI')
        do iprg = 1, prgtab%count()
            call prgtab%get_ui_program(prgnames(iprg), ptr2prg)
            call ptr2prg%create_json_entry(json, program_entry, ptr2prg%name%to_char(), 'program')
            call json%add(all_programs, program_entry)
        end do
        ! write & clean
        call json%print(all_programs, logfhandle)
        if ( json%failed() ) then
            THROW_HARD('json input/output error for simple_ui')
        endif
        call json%destroy(all_programs)
    end subroutine print_ui_json

    subroutine write_ui_json
        use json_module
        implicit none
        type(json_core)           :: json
        type(json_value), pointer :: all_programs, program_entry
        type(ui_program), pointer :: ptr2prg
        integer :: iprg
        ! JSON init
        call json%initialize(real_format=UI_JSON_REAL_FORMAT)
        ! create array of program entries
        call json%create_array(all_programs, 'SIMPLE User Interface')
        do iprg = 1, prgtab%count()
            call prgtab%get_ui_program(prgnames(iprg), ptr2prg)
            call ptr2prg%create_json_entry(json, program_entry, '', ptr2prg%name%to_char())
            call json%add(all_programs, program_entry)
        end do
        ! write & clean
        call json%print(all_programs, UI_FNAME)
        if ( json%failed() ) then
            THROW_HARD('json input/output error for simple_ui')
        endif
        call json%destroy(all_programs)
    end subroutine write_ui_json

    subroutine print_stream_ui_json
        ! this is ugly at the moment but paves the way ....
        use json_module
        type(json_core)           :: json
        type(json_value), pointer :: input, user_inputs, ui
        type(json_value), pointer :: processes, process, process_inputs
        ! JSON init
        call json%initialize(real_format=UI_JSON_REAL_FORMAT)
        ! create object of program entries
        call json%create_object(ui, 'STREAM_UI')
        ! master inputs
        call json%create_array(user_inputs, 'user_inputs')
        call json%add(ui, user_inputs)
        !! dir_movies
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'dir_movies')
        call json%add(input, 'keytype',     'dir')
        call json%add(input, 'label', 'Input movies directory')
        call json%add(input, 'help',  'Input movies directory')
        call json%add(input, 'required',    .TRUE.)
        !! dir_meta
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'dir_meta')
        call json%add(input, 'keytype',     'dir')
        call json%add(input, 'label', 'Input metadata directory')
        call json%add(input, 'help',  'Input metadata directory')
        call json%add(input, 'required',    .FALSE.)
        !! gainref
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'gainref')
        call json%add(input, 'keytype',     'file')
        call json%add(input, 'label', 'Gain reference')
        call json%add(input, 'help',  'Gain reference')
        call json%add(input, 'required',    .FALSE.)
        !! cs
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'cs')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'label', 'Spherical aberration (mm)')
        call json%add(input, 'help',  'Spherical aberration (mm)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     STREAM_DEFAULT_CS)
        !! fraca
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'fraca')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'label', 'Amplitude contrast fraction')
        call json%add(input, 'help',  'Amplitude contrast fraction')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     STREAM_DEFAULT_FRACA)
        !! kv
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'kv')
        call json%add(input, 'keytype',     'int')
        call json%add(input, 'label', 'Acceleration voltage (kV)')
        call json%add(input, 'help',  'Acceleration voltage (kV)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     int2str(STREAM_DEFAULT_KV))
        !! smpd
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'smpd')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'label', 'Pixel size (A)')
        call json%add(input, 'help',  'Pixel size (A)')
        call json%add(input, 'required',    .TRUE.)
        !! smpd_downscale
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'smpd_downscale')
        call json%add(input, 'keytype',     'hidden')
        call json%add(input, 'label', 'downscale pixel size (A)')
        call json%add(input, 'help',  'downscale pixel size (A)')
        call json%add(input, 'required',    .TRUE.)
        call json%add(input, 'default',     real2str(STREAM_DEFAULT_SMPD_DOWNSCALE))
        !! total_dose
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'total_dose')
        call json%add(input, 'keytype',     'float')
        call json%add(input, 'label', 'Total exposure dose (e/A2)')
        call json%add(input, 'help',  'Total exposure dose (e/A2)')
        call json%add(input, 'required',    .TRUE.)
        !! pickrefs
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'pickrefs')
        call json%add(input, 'keytype',     'file')
        call json%add(input, 'label', '2D averages for use as picking references (optional)')
        call json%add(input, 'help',  '2D averages for use as picking references (optional)')
        call json%add(input, 'required',    .FALSE.)
        !! box size
        call json%create_object(input, 'input')
        call json%add(user_inputs, input)
        call json%add(input, 'key',         'box_extract')
        call json%add(input, 'keytype',     'int')
        call json%add(input, 'label', 'force box size (px, optional)')
        call json%add(input, 'help',  'force a box size (px) eg. to match an existing dataset')
        call json%add(input, 'required',    .FALSE.)
        ! programs
        call json%create_array(processes, 'processes')
        call json%add(ui, processes)
        !! preproc
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         PREPROC_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'preproc')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'user_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_movies')
        call json%add(process_inputs, '', 'dir_meta')
        call json%add(process_inputs, '', 'gainref')
        call json%add(process_inputs, '', 'cs')
        call json%add(process_inputs, '', 'fraca')
        call json%add(process_inputs, '', 'kv')
        call json%add(process_inputs, '', 'smpd')
        call json%add(process_inputs, '', 'scale')
        call json%add(process_inputs, '', 'total_dose')
        call json%add(process_inputs, '', 'smpd_downscale')
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'outdir='   // PREPROC_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nparts='   // int2str(PREPROC_NPARTS))
        call json%add(process_inputs, '', 'nthr='     // int2str(PREPROC_NTHR))
        call json%add(process_inputs, '', 'ninipick=' // int2str(PREPROC_NINIPICK))
        !! assign_optics
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         OPTICS_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'assign_optics')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target=' // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='     // OPTICS_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nthr='       // int2str(DEFAULT_NTHR_MASTER))
        !! opening 2D
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         OPENING2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'gen_pickrefs')
        call json%add(process, 'nthr_master',  OPENING2D_NTHR)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // OPENING2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nthr='          // int2str(OPENING2D_NTHR))
        !! reference_based_picking
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         REFPICK_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'pick_extract')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'user_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'pickrefs')
        call json%add(process_inputs, '', 'box_extract')
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // PREPROC_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // REFPICK_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'nparts='        // int2str(REFPICK_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(REFPICK_NTHR))
        !! particle_sieving
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         SIEVING_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'sieve_cavgs')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // REFPICK_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // SIEVING_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'ncls='          // int2str(SIEVING_NCLS))
        call json%add(process_inputs, '', 'nptcls_per_cls='// int2str(SIEVING_NPTCLS_PER_CLASS))
        call json%add(process_inputs, '', 'nchunks='       // int2str(SIEVING_NCHUNKS))
        call json%add(process_inputs, '', 'nparts='        // int2str(SIEVING_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(SIEVING_NTHR))
        call json%add(process_inputs, '', 'interactive=yes')
        !! 2D classification
        call json%create_object(process, 'process')
        call json%add(processes, process)
        call json%add(process, 'name',         CLASS2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process, 'prg',          'abinitio2D_stream')
        call json%add(process, 'nthr_master',  DEFAULT_NTHR_MASTER)
        call json%create_array(process_inputs, 'static_inputs')
        call json%add(process, process_inputs)
        call json%add(process_inputs, '', 'dir_target='    // SIEVING_JOB_NAME)
        call json%add(process_inputs, '', 'optics_dir=../' // OPTICS_JOB_NAME)
        call json%add(process_inputs, '', 'outdir='        // CLASS2D_JOB_NAME) !important - directory names and name must match between processes
        call json%add(process_inputs, '', 'ncls='          // int2str(CLASS2D_NCLS))
        call json%add(process_inputs, '', 'nparts='        // int2str(CLASS2D_NPARTS))
        call json%add(process_inputs, '', 'nthr='          // int2str(CLASS2D_NTHR))
        ! print & clean
        call json%print(ui, logfhandle)
        if( json%failed() )then
            THROW_HARD('json input/output error for simple_ui')
        endif
        call json%destroy(ui)
    end subroutine print_stream_ui_json

    subroutine validate_ui_json
        use json_module
        use json_kinds
        use json_file_module
        type(json_core)                       :: json
        type(json_value),             pointer :: p
        character(kind=CK,len=:), allocatable :: fname
        ! Builds & writes UI
        call make_ui
        write(*,*)'Constructed UI'
        call write_ui_json
        write(*,*)'Wrote UI'
        ! Parses all values, check for duplicate and invalid types
        fname = UI_FNAME
        call json%parse(fname, p)
        write(*,*)'Completed json_parse_file'
        ! Cleanup
        call json%destroy()
        nullify(p)
    end subroutine validate_ui_json

end module simple_ui
