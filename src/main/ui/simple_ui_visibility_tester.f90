!@descr: unit tests for the UI visibility levels, descriptors and registered program contracts (simple_ui_visibility, simple_ui_program, simple_ui)
! Visibility levels and names; parameter descriptors (labels, help, unit inference from the
! descriptor text, unit-aware and generated placeholders, structured choices, CLI summaries);
! program descriptors (display names, categories, input visibility, groups, activation,
! choice and placeholder overrides, requirement groups); the registered programs' categories,
! requirements and CLI summaries; the registered test programs; and the phase-shift contract
! of the five CTF-fitting programs (fit_phshift binary, default no; phshift_min/max/step).
module simple_ui_visibility_tester
use simple_test_utils
use simple_linked_list,   only: linked_list, list_iterator
use simple_string,        only: string
use simple_ui,            only: make_ui, make_test_ui, get_prg_ptr, get_test_prg_ptr, count_prgs_in_category
use simple_ui_param,      only: UI_PLACEHOLDER_MAX_LEN, ui_param
use simple_ui_program,    only: UI_DISPLAY_NAME_MAX_LEN, UI_FILE, UI_PARM, UI_SUMMARY_MAX_LEN, &
    &category_descriptor, ui_cli_param_choices, ui_cli_param_summary, ui_program, ui_program_input, &
    &ui_activation_equals_any
use simple_ui_descriptor_types, only: ui_choices
use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_ADVANCED, UI_VIS_DEVELOPER, &
                               &ui_visibility_is_valid, ui_visibility_name
implicit none
private
public :: run_all_ui_visibility_tests

! shared by the lookup helpers below
type(ui_program), pointer :: registered_prg => null()
type(string) :: cli_text, program_name

contains

    subroutine run_all_ui_visibility_tests()
        write(*,'(A)') '**** running all UI visibility tests ****'
        call test_visibility_levels()
        call test_param_descriptors()
        call test_program_descriptor()
        call test_registered_programs()
        call test_registered_test_programs()
        call test_phshift_contract()
    end subroutine run_all_ui_visibility_tests

    subroutine test_visibility_levels()
        write(*,'(A)') 'test_visibility_levels'
        call assert_true_all_valid
        call assert_char('standard',  ui_visibility_name(UI_VIS_STANDARD),  'standard visibility name')
        call assert_char('advanced',  ui_visibility_name(UI_VIS_ADVANCED),  'advanced visibility name')
        call assert_char('developer', ui_visibility_name(UI_VIS_DEVELOPER), 'developer visibility name')
    end subroutine test_visibility_levels

    !> labels, help, unit inference, placeholders, choices and CLI summaries of single parameters
    subroutine test_param_descriptors()
        type(ui_param) :: param
        write(*,'(A)') 'test_param_descriptors'
        call param%set_param('test_param', 'num', 'Test parameter', 'Test parameter help', 'e.g. 1', .false., 1.)
        call assert_char('Test parameter', param%label%to_char(), 'input short text is the label')
        call assert_char('Test parameter help', param%help%to_char(), 'input long text is help')
        call assert_char('e.g. 10', param%placeholder%to_char(), 'numeric placeholder is a standard example')
        call assert_true(len_trim(param%placeholder%to_char()) <= UI_PLACEHOLDER_MAX_LEN, &
            &'placeholder respects the display limit')
        call assert_true(param%has_default, 'optional numeric parameter has a default')
        call assert_int(0, size(param%choices), 'numeric parameter has no choices')

        call param%set_param('mskdiam', 'num', 'Mask diameter', 'Mask diameter (in Angstroms)', 'e.g. 1', .true., 0.)
        call assert_char('Angstroms', param%units%to_char(), 'mask diameter unit is inferred from descriptor text')
        call assert_char('e.g. 180 Angstroms', param%placeholder%to_char(), 'mask diameter has a unit-aware example')
        cli_text = ui_cli_param_summary(param)
        call assert_char('Mask diameter; e.g. 180 Angstroms', cli_text%to_char(), 'CLI numeric summary includes the unit-aware example')

        call param%set_param('winsz', 'num', 'Window size', 'Window size (in pixels)', 'e.g. 1', .false., 1.)
        call assert_char('pixels', param%units%to_char(), 'pixel unit is inferred from descriptor text')
        call assert_char('e.g. 10 pixels', param%placeholder%to_char(), 'pixel placeholder includes its unit')
        call param%set_param('walltime', 'num', 'Walltime', 'Maximum execution time in seconds', 'e.g. 1', .false., 1.)
        call assert_char('seconds', param%units%to_char(), 'time unit is inferred from descriptor text')
        call assert_char('e.g. 10 seconds', param%placeholder%to_char(), 'time placeholder includes its unit')
        call param%set_param('bfac', 'num', 'B-factor', 'B-factor in Angstroms^2', 'e.g. 1', .false., 1.)
        call assert_char('Angstroms^2', param%units%to_char(), 'squared Angstrom unit is inferred from descriptor text')
        call assert_char('e.g. 10 Angstroms^2', param%placeholder%to_char(), 'squared Angstrom placeholder includes its unit')
        call param%set_param('cs', 'num', 'Spherical aberration', 'Spherical aberration constant (in mm)', 'e.g. 1', .true., 0.)
        call param%set_generated_default('2.7', .true.)
        call assert_char('mm', param%units%to_char(), 'millimeter unit is inferred from descriptor text')
        call assert_char('e.g. 2.7 mm', param%placeholder%to_char(), 'generated numeric default becomes the CLI example')
        call param%set_param('kv', 'num', 'Acceleration voltage', 'Acceleration voltage in kV', 'e.g. 1', .true., 0.)
        call param%set_generated_default('300.', .true.)
        call assert_char('e.g. 300 kV', param%placeholder%to_char(), 'generated numeric example is normalized and unit-aware')
        call param%set_param('nptcls_per_cls', 'num', 'Number of particles per cluster', 'Integer particle count', 'e.g. 1', .true., 0.)
        call param%set_generated_default('500', .true.)
        call assert_char('e.g. 500', param%placeholder%to_char(), 'integer generated default remains an integer CLI example')

        call param%set_param('subprojname', 'str', 'Subproject name', 'SIMPLE subproject name', 'e.g. value', .true., '')
        call assert_char('e.g. myproject.simple', param%placeholder%to_char(), &
            &'subproject names advertise the SIMPLE project extension')

        call param%set_param('projfile', 'file', 'Project file', 'SIMPLE project file', 'e.g. any-file', .true., '')
        call assert_char('e.g. input.simple', param%placeholder%to_char(), &
            &'project-file placeholder advertises the SIMPLE project format')

        call param%set_param('plaintexttab', 'file', 'Plain-text parameter table', 'Text parameters', 'e.g. any-file', .true., '')
        call assert_char('e.g. params.txt', param%placeholder%to_char(), &
            &'plain-text table placeholder advertises a text file')

        call param%set_param('stktab', 'file', 'Stack table', 'List of image stacks', 'e.g. any-file', .true., '')
        call assert_char('e.g. stktab.txt', param%placeholder%to_char(), &
            &'stack-table placeholder advertises a text file')

        call param%set_param('choice_param', 'multi', 'Choice parameter', 'Choice parameter help', &
            &'', .false., 'second', choices=ui_choices([character(len=6) :: 'first', 'second', 'third']))
        call assert_true(param%has_default, 'optional multi parameter has a default')
        call assert_int(3, size(param%choices), 'multi parameter choices are structured')
        call assert_char('first',  param%choices(1)%value%to_char(), 'first choice value')
        call assert_char('second', param%choices(2)%label%to_char(), 'choice label defaults to CLI value')
        call assert_char('', param%placeholder%to_char(), 'choice fields have no placeholder')
        cli_text = ui_cli_param_choices(param)
        call assert_char('(first|second|third)', cli_text%to_char(), 'CLI choice text uses structured values')
        cli_text = ui_cli_param_summary(param)
        call assert_char('Choice parameter (first|second|third){second}', cli_text%to_char(), &
            &'CLI choice summary includes the optional default')

        call param%set_param('binary_param', 'binary', 'Binary parameter', 'Binary parameter help', &
            &'', .true., '', choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call assert_true(.not. param%has_default, 'required binary parameter has no default')
        call assert_int(2, size(param%choices), 'binary parameter choices are structured')
        cli_text = ui_cli_param_summary(param)
        call assert_char('Binary parameter (yes|no)', cli_text%to_char(), 'CLI binary summary has no default marker')
    end subroutine test_param_descriptors

    !> a program descriptor: display name, category, input visibility, groups, activation, overrides, requirements
    subroutine test_program_descriptor()
        type(ui_param)   :: param
        type(ui_program) :: ui_prg
        type(string), allocatable :: supplied_keys(:)
        write(*,'(A)') 'test_program_descriptor'
        call param%set_param('binary_param', 'binary', 'Binary parameter', 'Binary parameter help', &
            &'', .true., '', choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
            &visibility=UI_VIS_STANDARD)
        call assert_char('Test program', ui_prg%summary%to_char(), 'program short text is the summary')
        call assert_true(len_trim(ui_prg%summary%to_char()) <= UI_SUMMARY_MAX_LEN, &
            &'program summary respects the display limit')
        call assert_char('Test program', ui_prg%display_name%to_char(), 'program display name falls back to the summary')
        call assert_true(len_trim(ui_prg%display_name%to_char()) <= UI_DISPLAY_NAME_MAX_LEN, &
            &'program display name respects the display limit')
        call assert_char('Test program help', ui_prg%help%to_char(), 'program long text is help')
        call assert_int(UI_VIS_STANDARD, ui_prg%visibility, 'explicit program standard visibility')
        call ui_prg%set_category(category_descriptor('test', 'Test Programs', 1))
        call assert_char('test', ui_prg%category%to_char(), 'program category')
        call assert_char('Test Programs', ui_prg%category_display_name%to_char(), 'program category display name')
        call assert_int(1, ui_prg%category_order, 'program category order')
        call ui_prg%new('test_program', 'Test program', 'Test program help', 'simple_exec', .false., &
            &visibility=UI_VIS_DEVELOPER, display_name='Visible test program')
        call assert_char('Visible test program', ui_prg%display_name%to_char(), 'explicit program display name')
        call assert_int(UI_VIS_DEVELOPER, ui_prg%visibility, 'explicit program developer visibility')
        call ui_prg%add_input(UI_PARM, 'visibility_param', 'num', 'Visibility parameter', &
            &'Visibility parameter help', 'e.g. 1', .false., 1., group='test controls', &
            &visibility=UI_VIS_STANDARD, activation=ui_activation_equals_any('mode', &
            &[character(len=8) :: 'standard', 'expert']))
        call assert_input_visibility(ui_prg%parm_ios, 'visibility_param', UI_VIS_STANDARD)
        call assert_input_binding(ui_prg%parm_ios, 'visibility_param')
        call ui_prg%add_input(UI_PARM, 'developer_param', 'num', 'Developer parameter', &
            &'Developer parameter help', 'e.g. 1', .false., 1., visibility=UI_VIS_DEVELOPER)
        call assert_input_visibility(ui_prg%parm_ios, 'developer_param', UI_VIS_DEVELOPER)
        call ui_prg%add_input(UI_PARM, param, choices_override=ui_choices([character(len=3) :: 'on', 'off']))
        call assert_input_visibility(ui_prg%parm_ios, 'binary_param', UI_VIS_STANDARD)
        call assert_input_choice(ui_prg%parm_ios, 'binary_param', 'off')
        call ui_prg%add_input(UI_FILE, 'input_a', 'file', 'Input A', 'First input source', 'e.g. input-a.mrc', .false., '')
        call ui_prg%add_input(UI_FILE, 'input_b', 'file', 'Input B', 'Second input source', 'e.g. input-b.mrc', .false., '')
        call assert_input_visibility(ui_prg%file_ios, 'input_a', UI_VIS_ADVANCED)
        call assert_input_placeholder(ui_prg%file_ios, 'input_a', 'e.g. input-a.mrc')
        call param%set_param('table', 'file', 'Input table', 'Input table help', 'e.g. table.txt', .false., '')
        call ui_prg%add_input(UI_FILE, param, placeholder_override='e.g. context.csv')
        call assert_input_placeholder(ui_prg%file_ios, 'table', 'e.g. context.csv')
        call ui_prg%add_requirement('input_source', 'Input source', 'Supply exactly one input source.', &
            &[character(len=7) :: 'input_a', 'input_b'], max_selected=1)
        call assert_int(1, size(ui_prg%requirements), 'one requirement group is registered')
        call assert_int(1, ui_prg%requirements(1)%min_selected, 'requirement minimum selection count')
        call assert_int(1, ui_prg%requirements(1)%max_selected, 'requirement maximum selection count')
        allocate(supplied_keys(1))
        supplied_keys(1) = 'input_a'
        call assert_true(ui_prg%requirements_satisfied(supplied_keys), 'one accepted input satisfies requirement')
        supplied_keys(1) = 'other'
        call assert_true(.not. ui_prg%requirements_satisfied(supplied_keys), 'unrelated input does not satisfy requirement')
        deallocate(supplied_keys)
    end subroutine test_program_descriptor

    !> categories, requirements and CLI summaries of the registered programs
    subroutine test_registered_programs()
        write(*,'(A)') 'test_registered_programs'
        call make_ui
        call assert_registered_requirement('binarize', 'input', 1, 1)
        call assert_registered_category('icm2D', 'denoise', 'Denoising', 70)
        call assert_registered_category('refine3D', 'refine3d', 'Refine 3D Workflows', 60)
        call assert_registered_category('flex_pca', 'heterogeneity', 'Heterogeneity Analysis', 65)
        call assert_registered_category('refine3D_states', 'heterogeneity', 'Heterogeneity Analysis', 65)
        call assert_registered_category('classify3D_refs', 'heterogeneity', 'Heterogeneity Analysis', 65)
        call assert_registered_category('ptcl3D_state_consensus', 'heterogeneity', 'Heterogeneity Analysis', 65)
        call assert_int(4, count_prgs_in_category('heterogeneity'), 'heterogeneity program count')
        call assert_registered_category('reconstruct3D', 'reconstruct3d', 'Reconstruct 3D Workflows', 68)
        call assert_registered_category('bootstrap_rec3D', 'reconstruct3d', 'Reconstruct 3D Workflows', 68)
        call assert_int(2, count_prgs_in_category('reconstruct3d'), 'reconstruct3d program count')
        call assert_registered_category('postprocess', 'postprocess', 'Post-processing', 69)
        call assert_registered_category('postprocess_nu', 'postprocess', 'Post-processing', 69)
        call assert_int(2, count_prgs_in_category('postprocess'), 'postprocess program count')
        call assert_registered_category('automask', 'mask', 'Masking', 100)
        call assert_registered_category('cls_split', 'cluster2d', 'Cluster2D Workflows', 30)
        call assert_registered_category('reimport_particles', 'project', 'Project Management', 10)
        call assert_registered_category('fractionate_movies', 'preproc', 'Pre-processing', 20)
        call assert_registered_category('split', 'image', 'General Image Processing', 90)
        call assert_registered_category('split_stack', 'image', 'General Image Processing', 90)
        call assert_registered_category('filter', 'filter', 'Filtering', 80)
        call assert_registered_category('new_project', 'project', 'Project Management', 10)
        call assert_registered_category('export_starproject', 'project', 'Project Management', 10)
        call assert_registered_category('stack', 'image', 'General Image Processing', 90)
        call assert_registered_category('motion_correct', 'preproc', 'Pre-processing', 20)
        call assert_program_not_registered('ppca_volvar')
        call assert_program_not_registered('export_manifoldem_starproject')
        call assert_registered_category('atoms_stats', 'atom', 'Atom Analysis', 50)
        call assert_registered_category('abinitio2D_stream', 'stream', 'Stream Workflows', 10)
        call assert_registered_category('abinitio2D', 'cluster2d', 'Cluster2D Workflows', 30)
        call assert_registered_cli_summary('model_cavgs_rejection', 'quality_mode', &
            &'Class-average quality mode (apply|analyze|learn|evaluate|promote){apply}')
        call assert_registered_cli_summary('import_movies', 'cs', 'Spherical aberration; e.g. 2.7 mm')
        call assert_registered_cli_summary('import_movies', 'fraca', 'Amplitude contrast fraction; e.g. 0.1')
        call assert_registered_cli_summary('import_movies', 'kv', 'Acceleration voltage; e.g. 300 kV')
        call assert_registered_cli_summary('extract', 'wfloat16', 'Write float16 particle stacks (yes|no){no}')
        call assert_registered_cli_summary('reextract', 'wfloat16', 'Write float16 particle stacks (yes|no){no}')
        call assert_registered_cli_summary('extract_subproj', 'subprojname', 'Subproject name; e.g. myproject.simple')
        call assert_registered_search_cli_summary('ctf_estimate', 'dfmin', 'Expected minimum defocus; e.g. 0.2 microns')
        call assert_registered_search_cli_summary('ctf_estimate', 'dfmax', 'Expected maximum defocus; e.g. 5.0 microns')
    end subroutine test_registered_programs

    !> the registered test programs
    subroutine test_registered_test_programs()
        write(*,'(A)') 'test_registered_test_programs'
        call make_test_ui
        program_name = 'cavg_registration'
        call get_test_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), 'cavg_registration test is registered')
        if( associated(registered_prg) )then
            call assert_char('utils', registered_prg%category%to_char(), 'cavg_registration test category')
        endif
        call assert_registered_test_category('preproc',    'stream', 'Stream',     130)
        call assert_registered_test_category('lib_stream', 'class',  'Unit tests', 10)
    end subroutine test_registered_test_programs

    !> the five CTF-fitting programs expose fit_phshift (binary, default no) with the phase-shift window
    subroutine test_phshift_contract()
        character(len=16), parameter :: FITTING_PROGRAMS(5) = [character(len=16) :: &
            &'ctf_estimate', 'preprocess', 'preproc', 'mini_stream', 'check_refpick']
        integer :: i
        write(*,'(A)') 'test_phshift_contract'
        call make_ui
        do i = 1, size(FITTING_PROGRAMS)
            program_name = trim(FITTING_PROGRAMS(i))
            call get_prg_ptr(program_name, registered_prg)
            call assert_true(associated(registered_prg), trim(FITTING_PROGRAMS(i))//' is registered')
            if( .not. associated(registered_prg) ) cycle
            call assert_ui_param(registered_prg%parm_ios, 'fit_phshift', FITTING_PROGRAMS(i), &
                &expected_type='binary', expected_default='no')
            call assert_ui_param(registered_prg%srch_ctrls, 'phshift_min',  FITTING_PROGRAMS(i))
            call assert_ui_param(registered_prg%srch_ctrls, 'phshift_max',  FITTING_PROGRAMS(i))
            call assert_ui_param(registered_prg%srch_ctrls, 'phshift_step', FITTING_PROGRAMS(i))
        enddo
    end subroutine test_phshift_contract

    subroutine assert_true_all_valid
        call assert_true(ui_visibility_is_valid(UI_VIS_STANDARD),  'standard visibility is valid')
        call assert_true(ui_visibility_is_valid(UI_VIS_ADVANCED),  'advanced visibility is valid')
        call assert_true(ui_visibility_is_valid(UI_VIS_DEVELOPER), 'developer visibility is valid')
        call assert_true(.not. ui_visibility_is_valid(0),          'zero visibility is invalid')
    end subroutine assert_true_all_valid

    subroutine assert_registered_category( name, expected_category, expected_display_name, expected_order )
        character(len=*), intent(in) :: name, expected_category, expected_display_name
        integer,          intent(in) :: expected_order
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if( associated(registered_prg) )then
            call assert_char(expected_category, registered_prg%category%to_char(), trim(name)//' category')
            call assert_char(expected_display_name, registered_prg%category_display_name%to_char(), &
                &trim(name)//' category display name')
            call assert_int(expected_order, registered_prg%category_order, trim(name)//' category order')
        endif
    end subroutine assert_registered_category

    subroutine assert_program_not_registered( name )
        character(len=*), intent(in) :: name
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(.not. associated(registered_prg), trim(name)//' is not registered')
    end subroutine assert_program_not_registered

    subroutine assert_registered_test_category( name, expected_category, expected_display_name, expected_order )
        character(len=*), intent(in) :: name, expected_category, expected_display_name
        integer,          intent(in) :: expected_order
        program_name = name
        call get_test_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' test is registered')
        if( associated(registered_prg) )then
            call assert_char(expected_category, registered_prg%category%to_char(), trim(name)//' test category')
            call assert_char(expected_display_name, registered_prg%category_display_name%to_char(), &
                &trim(name)//' test category display name')
            call assert_int(expected_order, registered_prg%category_order, trim(name)//' test category order')
            call assert_char('simple_test_exec', registered_prg%executable%to_char(), trim(name)//' test executable')
            call assert_int(0, registered_prg%get_nrequired_keys(), trim(name)//' test has no required inputs')
        endif
    end subroutine assert_registered_test_category

    subroutine assert_registered_requirement(name, expected_id, expected_minimum, expected_maximum)
        character(len=*), intent(in) :: name, expected_id
        integer, intent(in) :: expected_minimum, expected_maximum
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if (associated(registered_prg)) then
            call assert_int(1, size(registered_prg%requirements), trim(name)//' requirement count')
            if (size(registered_prg%requirements) == 1) then
                call assert_char(expected_id, registered_prg%requirements(1)%id%to_char(), &
                    &trim(name)//' requirement id')
                call assert_int(expected_minimum, registered_prg%requirements(1)%min_selected, &
                    &trim(name)//' requirement minimum')
                call assert_int(expected_maximum, registered_prg%requirements(1)%max_selected, &
                    &trim(name)//' requirement maximum')
            endif
        endif
    end subroutine assert_registered_requirement

    subroutine assert_registered_cli_summary(name, key, expected_summary)
        character(len=*), intent(in) :: name, key, expected_summary
        type(list_iterator)   :: iterator
        class(*), allocatable :: value
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if( .not. associated(registered_prg) ) return
        iterator = registered_prg%parm_ios%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        cli_text = ui_cli_param_summary(input%param)
                        call assert_char(expected_summary, cli_text%to_char(), trim(name)//' CLI choice summary')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, trim(name)//' input exists')
    end subroutine assert_registered_cli_summary

    subroutine assert_registered_search_cli_summary(name, key, expected_summary)
        character(len=*), intent(in) :: name, key, expected_summary
        type(list_iterator)   :: iterator
        class(*), allocatable :: value
        program_name = name
        call get_prg_ptr(program_name, registered_prg)
        call assert_true(associated(registered_prg), trim(name)//' is registered')
        if( .not. associated(registered_prg) ) return
        iterator = registered_prg%srch_ctrls%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        cli_text = ui_cli_param_summary(input%param)
                        call assert_char(expected_summary, cli_text%to_char(), trim(name)//' CLI search summary')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, trim(name)//' search input exists')
    end subroutine assert_registered_search_cli_summary

    subroutine assert_input_visibility( inputs, key, expected_visibility )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key
        integer,           intent(in) :: expected_visibility
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_int(expected_visibility, input%visibility, 'explicit input visibility')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_program_input')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input exists')
    end subroutine assert_input_visibility

    subroutine assert_input_choice( inputs, key, expected_choice )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key, expected_choice
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_char(expected_choice, input%param%choices(2)%value%to_char(), &
                            &'explicit choice override replaces structured choices')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_program_input')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input exists')
    end subroutine assert_input_choice

    subroutine assert_input_placeholder( inputs, key, expected_placeholder )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key, expected_placeholder
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_char(expected_placeholder, input%param%placeholder%to_char(), &
                            &'contextual input placeholder')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
                class default
                    call assert_int(1, 0, 'UI input is a ui_program_input')
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input placeholder exists')
    end subroutine assert_input_placeholder

    subroutine assert_input_binding( inputs, key )
        type(linked_list), intent(in) :: inputs
        character(len=*),  intent(in) :: key
        type(list_iterator)   :: iterator
        class(*), allocatable :: value

        iterator = inputs%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( input => value )
                type is( ui_program_input )
                    if( input%param%key%to_char() == key )then
                        call assert_char('test controls', input%group%id%to_char(), 'input group id')
                        call assert_int(1, input%group%order, 'first input group order')
                        call assert_char('mode', input%activation%key%to_char(), 'activation key')
                        call assert_int(2, size(input%activation%equals_any), 'activation value count')
                        call assert_char('expert', input%activation%equals_any(2)%to_char(), 'activation value')
                        if( allocated(value) ) deallocate(value)
                        return
                    endif
            end select
            if( allocated(value) ) deallocate(value)
            call iterator%next()
        enddo
        call assert_int(1, 0, 'UI input binding exists')
    end subroutine assert_input_binding

    !> the input `key` exists in `params`, with the expected type and default when given
    subroutine assert_ui_param( params, key, prg_name, expected_type, expected_default )
        type(linked_list), intent(in) :: params
        character(len=*),  intent(in) :: key, prg_name
        character(len=*),  intent(in), optional :: expected_type, expected_default
        type(list_iterator)   :: iterator
        class(*), allocatable :: value
        logical :: found
        found    = .false.
        iterator = params%begin()
        do while( iterator%has_value() )
            call iterator%getter(value)
            select type( param => value )
                type is( ui_program_input )
                    if( param%param%key%to_char() == key )then
                        found = .true.
                        if( present(expected_type) )then
                            call assert_char(expected_type, param%param%keytype%to_char(), &
                                &trim(prg_name)//': '//key//' UI type')
                        endif
                        if( present(expected_default) )then
                            call assert_char(expected_default, param%param%cval_default%to_char(), &
                                &trim(prg_name)//': '//key//' default')
                        endif
                    endif
                class default
                    call assert_true(.false., trim(prg_name)//': UI parameter-list entry is a ui_program_input')
            end select
            if( allocated(value) ) deallocate(value)
            if( found ) exit
            call iterator%next()
        enddo
        call assert_true(found, trim(prg_name)//': UI parameter '//key//' exists')
    end subroutine assert_ui_param

end module simple_ui_visibility_tester
