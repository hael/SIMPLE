!@descr: module defining the main user interface program type and associated methods
module simple_ui_program
use simple_core_module_api
use simple_ansi_ctrls
use simple_linked_list, only: linked_list, list_iterator
use simple_ui_param,    only: ui_param
use simple_ui_descriptor_types, only: ui_choice
use simple_ui_visibility, only: UI_VIS_STANDARD, UI_VIS_ADVANCED, UI_VIS_DEVELOPER, &
                               &ui_visibility_is_valid, ui_visibility_name
use simple_ui_default_values, only: get_ui_default
implicit none
#include "simple_local_flags.inc"

integer, parameter :: UI_IMG=1, UI_PARM=2, UI_FILE=3, UI_SRCH=4, UI_FILT=5, UI_MASK=6, UI_COMP=7
integer, parameter :: UI_SUMMARY_MIN_LEN = 30
integer, parameter :: UI_SUMMARY_MAX_LEN = 100
integer, parameter :: UI_DISPLAY_NAME_MAX_LEN = 100
character(len=*), parameter :: UI_JSON_REAL_FORMAT = '(ss,G0.6)'

type :: category_descriptor
    character(len=32) :: id = ''
    character(len=64) :: display_name = ''
    integer           :: order = 0
end type category_descriptor

type :: ui_input_group
    type(string) :: id
    type(string) :: label
    integer      :: order = 0
end type ui_input_group

type :: ui_activation
    type(string) :: key
    type(string), allocatable :: equals_any(:)
end type ui_activation

type :: ui_requirement_group
    type(string) :: id
    type(string) :: label
    type(string) :: help
    type(string), allocatable :: keys(:)
    integer :: min_selected = 1
    integer :: max_selected = 1
end type ui_requirement_group

type :: ui_program_input
    type(ui_param)       :: param
    integer              :: section = 0
    type(ui_input_group) :: group
    integer              :: visibility = UI_VIS_DEVELOPER
    type(ui_activation)  :: activation
end type ui_program_input

! production-level program interface for simple_exec, single_exec & simple_stream executables
type :: ui_program
    type(string) :: name
    type(string) :: category
    type(string) :: category_display_name
    integer      :: category_order = 0
    type(string) :: display_name
    type(string) :: summary
    type(string) :: help
    type(string) :: executable
    integer      :: visibility = UI_VIS_DEVELOPER
    type(ui_input_group), allocatable :: groups(:)
    type(ui_requirement_group), allocatable :: requirements(:)
    ! image input/output
    type(linked_list) :: img_ios
    ! file and directory input/output
    type(linked_list) :: file_ios
    ! parameter input/output
    type(linked_list) :: parm_ios
    ! search controls
    type(linked_list) :: srch_ctrls
    ! filter controls
    type(linked_list) :: filt_ctrls
    ! mask controls
    type(linked_list) :: mask_ctrls
    ! computer controls
    type(linked_list) :: comp_ctrls
    ! sp_project required flag
    logical :: sp_required = .true.
    ! existence flag
    logical :: exists = .false.
  contains
    procedure          :: new
    procedure, private :: add_input_num
    procedure, private :: add_input_str
    procedure, private :: add_input_param
    generic            :: add_input => add_input_num, add_input_str, add_input_param
    procedure          :: print_ui
    procedure          :: print_cmdline
    procedure          :: print_help
    procedure          :: create_json_entry
    procedure          :: write2json
    procedure          :: get_name
    procedure          :: get_category
    procedure          :: get_display_name
    procedure          :: get_executable
    procedure          :: set_category
    procedure          :: add_requirement
    procedure          :: requirements_satisfied
    procedure          :: get_nrequired_keys
    procedure          :: get_required_keys
    procedure          :: requires_sp_project
    procedure, private :: kill
end type ui_program

contains

    subroutine new( self, name, summary, help, executable, sp_required, visibility, display_name )
        class(ui_program),          intent(inout) :: self
        character(len=*),           intent(in)    :: name, summary, help, executable
        logical,                    intent(in)    :: sp_required
        integer,          optional, intent(in)    :: visibility
        character(len=*), optional, intent(in)    :: display_name
        call self%kill
        self%name        = trim(name)
        if( len_trim(summary) > UI_SUMMARY_MAX_LEN )then
            THROW_HARD('ui_program%new summary exceeds 100 characters: '//trim(name))
        endif
        self%summary     = trim(summary)
        if( present(display_name) )then
            if( len_trim(display_name) == 0 )then
                THROW_HARD('ui_program%new display_name must not be empty: '//trim(name))
            endif
            if( len_trim(display_name) > UI_DISPLAY_NAME_MAX_LEN )then
                THROW_HARD('ui_program%new display_name exceeds 100 characters: '//trim(name))
            endif
            self%display_name = trim(display_name)
        else
            ! Temporary migration fallback: every existing plain-English summary is a safe GUI title.
            self%display_name = trim(summary)
        endif
        self%help        = trim(help)
        self%executable  = trim(executable)
        self%sp_required = sp_required
        self%exists      = .true.
        if( present(visibility) )then
            if( .not. ui_visibility_is_valid(visibility) )then
                THROW_HARD('ui_program%new received an invalid visibility level')
            endif
            self%visibility = visibility
        endif
    end subroutine new

    subroutine add_input_num( self, which, key, keytype, label, help, placeholder, required, default_value, &
                            group, visibility, activation, choices, preserve_default )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        character(len=*),           intent(in)    :: key, keytype, label, help, placeholder
        logical,                    intent(in)    :: required
        real,                       intent(in)    :: default_value
        character(len=*), optional, intent(in)    :: group
        integer,          optional, intent(in)    :: visibility
        type(ui_activation), optional, intent(in) :: activation
        type(ui_choice), optional, intent(in) :: choices(:)
        logical, optional, intent(in) :: preserve_default
        type(ui_param)         :: p
        type(ui_program_input) :: input
        call p%set_param(key, keytype, label, help, placeholder, required, default_value, choices)
        call apply_generated_default(self, p, preserve_default)
        call bind_input(self, input, which, p, group, visibility, activation)
        call append_input(self, input)
    end subroutine add_input_num

    subroutine add_input_str( self, which, key, keytype, label, help, placeholder, required, default_value, &
                            group, visibility, activation, choices, preserve_default )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        character(len=*),           intent(in)    :: key, keytype, label, help, placeholder
        logical,                    intent(in)    :: required
        character(len=*),           intent(in)    :: default_value
        character(len=*), optional, intent(in)    :: group
        integer,          optional, intent(in)    :: visibility
        type(ui_activation), optional, intent(in) :: activation
        type(ui_choice), optional, intent(in) :: choices(:)
        logical, optional, intent(in) :: preserve_default
        type(ui_param)         :: p
        type(ui_program_input) :: input
        call p%set_param(key, keytype, label, help, placeholder, required, default_value, choices)
        call apply_generated_default(self, p, preserve_default)
        call bind_input(self, input, which, p, group, visibility, activation)
        call append_input(self, input)
    end subroutine add_input_str

    subroutine add_input_param( self, which, param, label_override, help_override, placeholder_override,&
    &required_override, group, visibility, activation, choices_override, preserve_default )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        type(ui_param),             intent(in)    :: param
        character(len=*), optional, intent(in)    :: label_override, help_override, placeholder_override
        logical,          optional, intent(in)    :: required_override
        character(len=*), optional, intent(in)    :: group
        integer,          optional, intent(in)    :: visibility
        type(ui_activation), optional, intent(in) :: activation
        type(ui_choice), optional, intent(in) :: choices_override(:)
        logical, optional, intent(in) :: preserve_default
        type(ui_param)         :: p
        type(ui_program_input) :: input
        p = param
        if( present(label_override)       ) p%label       = label_override
        if( present(help_override)        ) p%help        = help_override
        if( present(choices_override) ) call p%set_choices(choices_override)
        if( present(placeholder_override) ) call p%set_placeholder(placeholder_override)
        if( present(required_override) )then
            p%required = required_override
            if( p%required ) p%has_default = .false.
        endif
        call apply_generated_default(self, p, preserve_default)
        call bind_input(self, input, which, p, group, visibility, activation)
        call append_input(self, input)
    end subroutine add_input_param

    subroutine apply_generated_default(self, param, preserve_default)
        class(ui_program), intent(in)    :: self
        type(ui_param),    intent(inout) :: param
        logical, optional,  intent(in)    :: preserve_default
        character(len=128) :: value
        logical            :: found
        if( present(preserve_default) )then
            if( preserve_default ) return
        endif
        call get_ui_default(self%executable%to_char(), self%name%to_char(), param%key%to_char(), found, value)
        call param%set_generated_default(value, found)
    end subroutine apply_generated_default

    function ui_activation_equals_any( key, values ) result( activation )
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: values(:)
        type(ui_activation) :: activation
        integer :: i
        if( len_trim(key) == 0 )then
            THROW_HARD('ui_activation_equals_any received an empty key')
        endif
        if( size(values) == 0 )then
            THROW_HARD('ui_activation_equals_any requires at least one value')
        endif
        activation%key = trim(key)
        allocate(activation%equals_any(size(values)))
        do i = 1, size(values)
            if( len_trim(values(i)) == 0 )then
                THROW_HARD('ui_activation_equals_any received an empty value')
            endif
            activation%equals_any(i) = trim(values(i))
        enddo
    end function ui_activation_equals_any

    subroutine bind_input( self, input, section, param, group, visibility, activation )
        class(ui_program),          intent(inout) :: self
        type(ui_program_input),     intent(out)   :: input
        integer,                    intent(in)    :: section
        type(ui_param),             intent(in)    :: param
        character(len=*), optional, intent(in)    :: group
        integer,          optional, intent(in)    :: visibility
        type(ui_activation), optional, intent(in) :: activation
        input%param   = param
        input%section = section
        if( present(visibility) )then
            if( .not. ui_visibility_is_valid(visibility) )then
                THROW_HARD('ui_program input received an invalid visibility level')
            endif
        endif
        ! Requiredness is final at this point and required inputs are always
        ! part of the Standard interface. Optional inputs are context-specific.
        if( input%param%required )then
            if( present(visibility) )then
                if( visibility /= UI_VIS_STANDARD )then
                    THROW_HARD('required ui_program input must have Standard visibility')
                endif
            endif
            input%visibility = UI_VIS_STANDARD
        else
            input%visibility = UI_VIS_ADVANCED
            if( present(visibility) ) input%visibility = visibility
        endif
        if( present(group) ) call register_group(self, group, input%group)
        if( present(activation) ) input%activation = activation
    end subroutine bind_input

    subroutine append_input( self, input )
        class(ui_program),      intent(inout) :: self
        type(ui_program_input), intent(in)    :: input
        select case (input%section)
            case (UI_IMG);  call self%img_ios%push_back(input)
            case (UI_FILE); call self%file_ios%push_back(input)
            case (UI_PARM); call self%parm_ios%push_back(input)
            case (UI_SRCH); call self%srch_ctrls%push_back(input)
            case (UI_FILT); call self%filt_ctrls%push_back(input)
            case (UI_MASK); call self%mask_ctrls%push_back(input)
            case (UI_COMP); call self%comp_ctrls%push_back(input)
            case default
                THROW_HARD('which field selector: '//int2str(input%section)//' is unsupported; ui_program::add_input')
        end select
    end subroutine append_input

    subroutine add_requirement(self, id, label, help, keys, min_selected, max_selected)
        class(ui_program), intent(inout) :: self
        character(len=*),  intent(in)    :: id, label, help
        character(len=*),  intent(in)    :: keys(:)
        integer, optional, intent(in)    :: min_selected, max_selected
        type(ui_requirement_group), allocatable :: requirements(:)
        integer :: i, j, min_keys, max_keys, nrequirements
        if (.not. self%exists) THROW_HARD('ui_program%add_requirement called before ui_program%new')
        if (len_trim(id) == 0 .or. len_trim(label) == 0 .or. len_trim(help) == 0) then
            THROW_HARD('ui_program%add_requirement requires id, label, and help')
        endif
        if (size(keys) == 0) THROW_HARD('ui_program%add_requirement requires at least one key')
        min_keys = 1
        if (present(min_selected)) min_keys = min_selected
        max_keys = size(keys)
        if (present(max_selected)) max_keys = max_selected
        if (min_keys < 1 .or. max_keys < min_keys .or. max_keys > size(keys)) then
            THROW_HARD('ui_program%add_requirement received invalid cardinality')
        endif
        do i = 1, size(keys)
            if (len_trim(keys(i)) == 0) THROW_HARD('ui_program%add_requirement received an empty key')
            if (.not. has_input_key(self, keys(i))) then
                THROW_HARD('ui_program%add_requirement key is not registered: '//trim(keys(i)))
            endif
            do j = 1, i - 1
                if (trim(keys(i)) == trim(keys(j))) THROW_HARD('ui_program%add_requirement has a duplicate key')
            enddo
        enddo
        nrequirements = 0
        if (allocated(self%requirements)) nrequirements = size(self%requirements)
        do i = 1, nrequirements
            if (self%requirements(i)%id == trim(id)) THROW_HARD('ui_program%add_requirement has a duplicate id')
        enddo
        allocate(requirements(nrequirements + 1))
        if (nrequirements > 0) requirements(:nrequirements) = self%requirements
        requirements(nrequirements + 1)%id           = trim(id)
        requirements(nrequirements + 1)%label        = trim(label)
        requirements(nrequirements + 1)%help         = trim(help)
        requirements(nrequirements + 1)%min_selected = min_keys
        requirements(nrequirements + 1)%max_selected = max_keys
        allocate(requirements(nrequirements + 1)%keys(size(keys)))
        do i = 1, size(keys)
            requirements(nrequirements + 1)%keys(i) = trim(keys(i))
        enddo
        call move_alloc(requirements, self%requirements)
    end subroutine add_requirement

    logical function requirements_satisfied(self, defined_keys)
        class(ui_program), intent(in) :: self
        type(string),      intent(in) :: defined_keys(:)
        integer :: i, j, nselected
        requirements_satisfied = .true.
        if (.not. allocated(self%requirements)) return
        do i = 1, size(self%requirements)
            nselected = 0
            do j = 1, size(self%requirements(i)%keys)
                if (key_is_defined(self%requirements(i)%keys(j), defined_keys)) nselected = nselected + 1
            enddo
            if (nselected < self%requirements(i)%min_selected .or. &
                nselected > self%requirements(i)%max_selected) then
                requirements_satisfied = .false.
                return
            endif
        enddo
    end function requirements_satisfied

    subroutine register_group( self, group_id, group )
        class(ui_program),      intent(inout) :: self
        character(len=*),       intent(in)    :: group_id
        type(ui_input_group),   intent(out)   :: group
        type(ui_input_group), allocatable :: expanded_groups(:)
        integer :: i, ngroups
        if( len_trim(group_id) == 0 )then
            THROW_HARD('ui_program%register_group received an empty group id')
        endif
        if( allocated(self%groups) )then
            do i = 1, size(self%groups)
                if( trim(self%groups(i)%id%to_char()) == trim(group_id) )then
                    group = self%groups(i)
                    return
                endif
            enddo
        endif
        ngroups = 1
        if( allocated(self%groups) ) ngroups = size(self%groups) + 1
        allocate(expanded_groups(ngroups))
        if( allocated(self%groups) ) expanded_groups(:ngroups - 1) = self%groups
        expanded_groups(ngroups)%id    = trim(group_id)
        expanded_groups(ngroups)%label = group_label(group_id)
        expanded_groups(ngroups)%order = ngroups
        call move_alloc(expanded_groups, self%groups)
        group = self%groups(ngroups)
    end subroutine register_group

    pure function group_label( group_id ) result( label )
        character(len=*), intent(in) :: group_id
        character(len=64) :: label
        select case( trim(group_id) )
            case('data');              label = 'Data'
            case('model');             label = 'Model'
            case('search');            label = 'Search'
            case('filter');            label = 'Filter'
            case('mask');              label = 'Mask'
            case('compute');           label = 'Compute'
            case('quality');           label = 'Quality'
            case('picking');           label = 'Picking'
            case('extract');           label = 'Extract'
            case('motion correction'); label = 'Motion Correction'
            case('CTF estimation');    label = 'CTF Estimation'
            case('optics groups');     label = 'Optics Groups'
            case('cluster 2D');        label = '2D Clustering'
            case default;               label = trim(group_id)
        end select
    end function group_label

    subroutine add_program_groups_json( json, program, groups )
        use json_module
        class(json_core),                    intent(inout) :: json
        type(json_value), pointer,           intent(inout) :: program
        type(ui_input_group), allocatable,   intent(in)    :: groups(:)
        type(json_value), pointer :: group_entries, group_entry
        integer :: i
        if( .not. allocated(groups) ) return
        call json%create_array(group_entries, 'groups')
        do i = 1, size(groups)
            call json%create_object(group_entry, '')
            call json%add(group_entry, 'id', groups(i)%id%to_char())
            call json%add(group_entry, 'label', groups(i)%label%to_char())
            call json%add(group_entry, 'order', groups(i)%order)
            call json%add(group_entries, group_entry)
        enddo
        call json%add(program, group_entries)
    end subroutine add_program_groups_json

    subroutine add_program_input_json( json, entry, input )
        use json_module
        class(json_core),          intent(inout) :: json
        type(json_value), pointer, intent(inout) :: entry
        type(ui_program_input),    intent(in)    :: input
        type(json_value), pointer :: group_entry, activation_entry, activation_values, options
        integer :: i
        call json%add(entry, 'key',         input%param%key%to_char())
        call json%add(entry, 'keytype',     input%param%keytype%to_char())
        call json%add(entry, 'label',       input%param%label%to_char())
        call json%add(entry, 'help',        input%param%help%to_char())
        call json%add(entry, 'placeholder', input%param%placeholder%to_char())
        call json%add(entry, 'required',    input%param%required)
        call json%add(entry, 'has_default', input%param%has_default)
        call json%add(entry, 'visibility', trim(ui_visibility_name(input%visibility)))
        if( len_trim(input%param%units%to_char()) > 0 )then
            call json%add(entry, 'units', input%param%units%to_char())
        endif
        if( input%param%has_default )then
            if( input%param%keytype%to_char() == 'num' .or. input%param%keytype%to_char() == 'int' .or. &
                &input%param%keytype%to_char() == 'float' )then
                call json%add(entry, 'default', real(input%param%rval_default,dp))
            else
                call json%add(entry, 'default', input%param%cval_default%to_char())
            endif
        endif
        if( input%group%order > 0 )then
            call json%create_object(group_entry, 'group')
            call json%add(group_entry, 'id', input%group%id%to_char())
            call json%add(group_entry, 'label', input%group%label%to_char())
            call json%add(group_entry, 'order', input%group%order)
            call json%add(entry, group_entry)
        endif
        if( allocated(input%activation%equals_any) )then
            call json%create_object(activation_entry, 'activation')
            call json%add(activation_entry, 'key', input%activation%key%to_char())
            call json%create_array(activation_values, 'equals_any')
            do i = 1, size(input%activation%equals_any)
                call json%add(activation_values, '', input%activation%equals_any(i)%to_char())
            enddo
            call json%add(activation_entry, activation_values)
            call json%add(entry, activation_entry)
        endif
        if( allocated(input%param%choices) )then
            call json%create_array(options, 'options')
            do i = 1, size(input%param%choices)
                call json%add(options, '', input%param%choices(i)%value%to_char())
            enddo
            call json%add(entry, options)
        endif
    end subroutine add_program_input_json

    subroutine print_ui( self )
        class(ui_program), intent(in) :: self
        type(chash) :: ch
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') '>>> PROGRAM INFO'
        call ch%new(9)
        call ch%push('name',         self%name%to_char())
        call ch%push('category',     self%category%to_char())
        call ch%push('category_name', self%category_display_name%to_char())
        call ch%push('category_order', int2str(self%category_order))
        call ch%push('display_name', self%display_name%to_char())
        call ch%push('summary',     self%summary%to_char())
        call ch%push('help',        self%help%to_char())
        call ch%push('executable',  self%executable%to_char())
        call ch%push('visibility',  trim(ui_visibility_name(self%visibility)))
        call ch%print_key_val_pairs(logfhandle)
        call ch%kill
        if (.not. self%img_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
            call print_param_list(self%img_ios)
        endif
        if (.not. self%file_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('FILE INPUT/OUTPUT',      C_UNDERLINED)
            call print_param_list(self%file_ios)
        endif
        if (.not. self%parm_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
            call print_param_list(self%parm_ios)
        endif
        if (.not. self%srch_ctrls%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('SEARCH CONTROLS',        C_UNDERLINED)
            call print_param_list(self%srch_ctrls)
        endif
        if (.not. self%filt_ctrls%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('FILTER CONTROLS',        C_UNDERLINED)
            call print_param_list(self%filt_ctrls)
        endif
        if (.not. self%mask_ctrls%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('MASK CONTROLS',          C_UNDERLINED)
            call print_param_list(self%mask_ctrls)
        endif
        if (.not. self%comp_ctrls%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('COMPUTER CONTROLS',      C_UNDERLINED)
            call print_param_list(self%comp_ctrls)
        endif
    end subroutine print_ui

    subroutine print_cmdline( self, defined_keys )
        class(ui_program), intent(in) :: self
        type(string), optional, intent(in) :: defined_keys(:)
        ! type(string), intent(in) :: exec_cmd
        write(logfhandle,'(a)') format_str('USAGE', C_UNDERLINED)
        ! select case(exec_cmd%to_char())
        !     case('simple_exec')
        !         write(logfhandle,'(a)') format_str('bash-3.2$ simple_exec prg=' // self%name%to_char() // ' key1=val1 key2=val2 ...', C_ITALIC)
        !     case('simple_test_exec')
        !         write(logfhandle,'(a)') format_str('bash-3.2$ simple_test_exec prg=' // self%name%to_char() // ' key1=val1 key2=val2 ...', C_ITALIC)
        !     case('single_exec')
        !         write(logfhandle,'(a)') format_str('bash-3.2$ single_exec prg=' // self%name%to_char() // ' key1=val1 key2=val2 ...', C_ITALIC)
        !     case DEFAULT
        !         write(logfhandle,'(a)') format_str('bash-3.2$ ' // exec_cmd%to_char() // ' prg=' // self%name%to_char() // ' key1=val1 key2=val2 ...', C_ITALIC)
        ! end select
        write(logfhandle,'(a)') 'Required input parameters are shown in ' // format_str('bold', C_BOLD) // &
            &'; additional input requirements are listed below.'
        if (.not. self%img_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT', C_UNDERLINED)
            call print_param_hash(self%img_ios)
        end if
        if (.not. self%file_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('FILE INPUT/OUTPUT', C_UNDERLINED)
            call print_param_hash(self%file_ios)
        end if
        if (.not. self%parm_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
            call print_param_hash(self%parm_ios)
        end if
        if (.not. self%srch_ctrls%is_empty()) then
            write(logfhandle,'(a)') format_str('SEARCH CONTROLS', C_UNDERLINED)
            call print_param_hash(self%srch_ctrls)
        end if
        if (.not. self%filt_ctrls%is_empty()) then
            write(logfhandle,'(a)') format_str('FILTER CONTROLS', C_UNDERLINED)
            call print_param_hash(self%filt_ctrls)
        end if
        if (.not. self%mask_ctrls%is_empty()) then
            write(logfhandle,'(a)') format_str('MASK CONTROLS', C_UNDERLINED)
            call print_param_hash(self%mask_ctrls)
        end if
        if (.not. self%comp_ctrls%is_empty()) then
            write(logfhandle,'(a)') format_str('COMPUTER CONTROLS', C_UNDERLINED)
            call print_param_hash(self%comp_ctrls)
        end if
        call print_requirement_groups(self, defined_keys)
    end subroutine print_cmdline

    subroutine print_help( self )
        class(ui_program), intent(in) :: self
        write(logfhandle,'(a)') self%help%to_char()
    end subroutine print_help

    subroutine write2json( self )
        use json_module
        class(ui_program), intent(in) :: self
        type(json_core)               :: json
        type(json_value), pointer     :: program_entry
        ! JSON init
        call json%initialize(real_format=UI_JSON_REAL_FORMAT)
        call self%create_json_entry(json, program_entry, '', self%name%to_char(), .false.)
        ! write & clean
        call json%print(program_entry, self%name%to_char()//'.json')
        if( json%failed() )then
            THROW_HARD('json input/output error for program: '//self%name%to_char())
        endif
        call json%destroy(program_entry)
    end subroutine write2json

    subroutine create_json_entry(self, json, program_entry, entry_name, program_name, requirements_first)
        use json_module
        class(ui_program), intent(in) :: self
        class(json_core), intent(inout) :: json
        type(json_value), pointer, intent(out) :: program_entry
        character(len=*), intent(in) :: entry_name, program_name
        logical, optional, intent(in) :: requirements_first
        type(json_value), pointer :: program
        logical :: add_requirements_first

        add_requirements_first = .true.
        if (present(requirements_first)) add_requirements_first = requirements_first
        call json%create_object(program_entry, entry_name)
        call json%create_object(program, program_name)
        call json%add(program_entry, program)
        call json%add(program, 'name', self%name%to_char())
        call json%add(program, 'category', self%category%to_char())
        call json%add(program, 'category_display_name', self%category_display_name%to_char())
        call json%add(program, 'category_order', self%category_order)
        call json%add(program, 'display_name', self%display_name%to_char())
        call json%add(program, 'summary', self%summary%to_char())
        call json%add(program, 'help', self%help%to_char())
        call json%add(program, 'executable', self%executable%to_char())
        call json%add(program, 'visibility', trim(ui_visibility_name(self%visibility)))
        call add_program_groups_json(json, program, self%groups)
        if (add_requirements_first) call add_program_requirements_json(json, program, self%requirements)
        call create_section_from_list(json, program_entry, 'image input/output', self%img_ios)
        call create_section_from_list(json, program_entry, 'file input/output', self%file_ios)
        call create_section_from_list(json, program_entry, 'parameter input/output', self%parm_ios)
        call create_section_from_list(json, program_entry, 'search controls', self%srch_ctrls)
        call create_section_from_list(json, program_entry, 'filter controls', self%filt_ctrls)
        call create_section_from_list(json, program_entry, 'mask controls', self%mask_ctrls)
        call create_section_from_list(json, program_entry, 'computer controls', self%comp_ctrls)
        if (.not. add_requirements_first) call add_program_requirements_json(json, program, self%requirements)
    end subroutine create_json_entry

    function get_name( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%name
    end function get_name

    function get_category( self ) result( category )
        class(ui_program), intent(in) :: self
        type(string) :: category
        category = self%category
    end function get_category

    function get_display_name( self ) result( display_name )
        class(ui_program), intent(in) :: self
        type(string) :: display_name
        display_name = self%display_name
    end function get_display_name

    function get_executable( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%executable
    end function get_executable

    subroutine set_category( self, category )
        class(ui_program), intent(inout) :: self
        type(category_descriptor), intent(in) :: category
        if( len_trim(category%id) == 0 )then
            THROW_HARD('ui_program%set_category received an empty category')
        endif
        if( len_trim(category%display_name) == 0 )then
            THROW_HARD('ui_program%set_category received an empty category display name')
        endif
        if( category%order <= 0 )then
            THROW_HARD('ui_program%set_category received a non-positive category order')
        endif
        self%category              = trim(category%id)
        self%category_display_name = trim(category%display_name)
        self%category_order        = category%order
    end subroutine set_category

    integer function get_nrequired_keys( self )
        class(ui_program), intent(in) :: self
        get_nrequired_keys = count_required_in_list(self%img_ios)   + &
                            count_required_in_list(self%file_ios)  + &
                            count_required_in_list(self%parm_ios)   + &
                            count_required_in_list(self%srch_ctrls) + &
                            count_required_in_list(self%filt_ctrls) + &
                            count_required_in_list(self%mask_ctrls) + &
                            count_required_in_list(self%comp_ctrls)
    end function get_nrequired_keys

    function get_required_keys( self ) result( keys )
        class(ui_program), intent(in) :: self
        type(string), allocatable     :: keys(:)
        integer                       :: nreq, ireq
        ! count # required
        nreq = self%get_nrequired_keys()
        if (nreq <= 0) return
        allocate(keys(nreq))
        ireq = 0
        call append_required_keys_from_list(self%img_ios,    keys, ireq)
        call append_required_keys_from_list(self%file_ios,   keys, ireq)
        call append_required_keys_from_list(self%parm_ios,   keys, ireq)
        call append_required_keys_from_list(self%srch_ctrls, keys, ireq)
        call append_required_keys_from_list(self%filt_ctrls, keys, ireq)
        call append_required_keys_from_list(self%mask_ctrls, keys, ireq)
        call append_required_keys_from_list(self%comp_ctrls, keys, ireq)
        ! shrink if needed
        if (ireq < nreq) keys = keys(:ireq)
    end function get_required_keys

    logical function requires_sp_project( self )
        class(ui_program), intent(in) :: self
        requires_sp_project = self%sp_required
    end function requires_sp_project

    subroutine kill( self )
        class(ui_program), intent(inout) :: self
        if (.not. self%exists) return
        call self%name%kill()
        call self%category%kill()
        call self%category_display_name%kill()
        call self%display_name%kill()
        call self%summary%kill()
        call self%help%kill()
        call self%executable%kill()
        if( allocated(self%groups) ) deallocate(self%groups)
        if( allocated(self%requirements) ) deallocate(self%requirements)
        call self%img_ios%kill()
        call self%file_ios%kill()
        call self%parm_ios%kill()
        call self%srch_ctrls%kill()
        call self%filt_ctrls%kill()
        call self%mask_ctrls%kill()
        call self%comp_ctrls%kill()
        self%visibility  = UI_VIS_DEVELOPER
        self%category_order = 0
        self%sp_required = .true.
        self%exists      = .false.
    end subroutine kill

    ! private helpers

    subroutine print_param_list( lst )
        class(linked_list), intent(in) :: lst
        type(list_iterator)   :: it
        class(*), allocatable :: tmp
        type(chash)           :: ch
        integer               :: i
        if (lst%is_empty()) return
        i  = 0
        it = lst%begin()
        do while (it%has_value())
            i = i + 1
            call it%getter(tmp)
            select type(t => tmp)
                type is (ui_program_input)
                    write(logfhandle,'(a,1x,i3)') '>>> PARAMETER #', i
                    call ch%new(8)
                    call ch%push('key',         t%param%key%to_char())
                    call ch%push('keytype',     t%param%keytype%to_char())
                    call ch%push('label',       t%param%label%to_char())
                    call ch%push('help',        t%param%help%to_char())
                    call ch%push('placeholder', t%param%placeholder%to_char())
                    call ch%push('required', merge('T','F', t%param%required))
                    call ch%push('visibility', trim(ui_visibility_name(t%visibility)))
                    call ch%push('group', t%group%label%to_char())
                    call ch%print_key_val_pairs(logfhandle)
                    call ch%kill
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('print_param_list: list element is not ui_program_input')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        end do
    end subroutine print_param_list

    subroutine print_param_hash( lst )
        class(linked_list),    intent(in)  :: lst
        character(len=KEYLEN), allocatable :: keys(:), sorted_keys(:), rearranged_keys(:)
        logical,               allocatable :: req(:), sorted_req(:)
        integer,               allocatable :: inds(:)
        type(chash)            :: ch
        type(list_iterator)    :: it
        class(*), allocatable  :: tmp
        integer :: i, nparams, nreq, iopt
        if (lst%is_empty()) return
        nparams = lst%size()
        call ch%new(nparams)
        allocate(keys(nparams), sorted_keys(nparams), rearranged_keys(nparams), req(nparams),  sorted_req(nparams))
        ! gather
        i  = 0
        it = lst%begin()
        do while (it%has_value())
            i = i + 1
            call it%getter(tmp)
            select type(t => tmp)
            type is (ui_program_input)
                call ch%push(t%param%key%to_char(), t%param%label%to_char()//'; '//t%param%placeholder%to_char())
                keys(i) = t%param%key%to_char()
                req(i)  = t%param%required
            class default
                if (allocated(tmp)) deallocate(tmp)
                THROW_HARD('print_param_hash: list element is not ui_program_input')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        end do
        ! sort keys, keep req aligned via inds
        sorted_keys = keys
        call lex_sort(sorted_keys, inds=inds)
        if (allocated(inds)) then
            do i = 1, nparams
                sorted_req(i) = req(inds(i))
            end do
        else
            ! if lex_sort doesn’t return inds for some reason, fall back
            sorted_req = req
        end if
        ! required-first reordering (within already-sorted order)
        if (any(sorted_req)) then
            nreq = 0
            do i = 1, nparams
                if (sorted_req(i)) then
                    nreq = nreq + 1
                    rearranged_keys(nreq) = sorted_keys(i)
                end if
            end do
            iopt = nreq
            do i = 1, nparams
                if (.not. sorted_req(i)) then
                    iopt = iopt + 1
                    rearranged_keys(iopt) = sorted_keys(i)
                end if
            end do
            sorted_keys = rearranged_keys
            sorted_req(:nreq)     = .true.
            sorted_req(nreq+1:)   = .false.
        end if
        call ch%print_key_val_pairs(logfhandle, sorted_keys, mask=sorted_req)
        call ch%kill()
        if (allocated(keys))            deallocate(keys)
        if (allocated(sorted_keys))     deallocate(sorted_keys)
        if (allocated(rearranged_keys)) deallocate(rearranged_keys)
        if (allocated(req))             deallocate(req)
        if (allocated(sorted_req))      deallocate(sorted_req)
        if (allocated(inds))            deallocate(inds)
    end subroutine print_param_hash

    subroutine create_section_from_list( json, program_entry, name, lst )
        use json_module
        class(json_core),          intent(inout) :: json
        type(json_value), pointer, intent(inout) :: program_entry
        character(len=*),          intent(in)    :: name
        class(linked_list),        intent(in)    :: lst
        type(json_value), pointer :: entry, section
        type(list_iterator)       :: it
        class(*), allocatable     :: tmp
        call json%create_array(section, trim(name))
        if (.not. lst%is_empty()) then
            it = lst%begin()
            do while (it%has_value())
                call it%getter(tmp)
                select type(t => tmp)
                type is (ui_program_input)
                    call json%create_object(entry, t%param%key%to_char())
                    call add_program_input_json(json, entry, t)
                    call json%add(section, entry)
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('create_section_from_list: list element is not ui_program_input')
                end select
                if (allocated(tmp)) deallocate(tmp)
                call it%next()
            end do
        end if
        call json%add(program_entry, section)
    end subroutine create_section_from_list

    subroutine add_program_requirements_json(json, program, requirements)
        use json_module
        class(json_core), intent(inout) :: json
        type(json_value), pointer, intent(inout) :: program
        type(ui_requirement_group), allocatable, intent(in) :: requirements(:)
        type(json_value), pointer :: requirement, requirement_list, keys
        integer :: i, j
        call json%create_array(requirement_list, 'requirements')
        if (allocated(requirements)) then
            do i = 1, size(requirements)
                call json%create_object(requirement, requirements(i)%id%to_char())
                call json%add(requirement, 'id', requirements(i)%id%to_char())
                call json%add(requirement, 'label', requirements(i)%label%to_char())
                call json%add(requirement, 'help', requirements(i)%help%to_char())
                call json%add(requirement, 'min_selected', requirements(i)%min_selected)
                call json%add(requirement, 'max_selected', requirements(i)%max_selected)
                call json%create_array(keys, 'keys')
                do j = 1, size(requirements(i)%keys)
                    call json%add(keys, '', requirements(i)%keys(j)%to_char())
                enddo
                call json%add(requirement, keys)
                call json%add(requirement_list, requirement)
            enddo
        endif
        call json%add(program, requirement_list)
    end subroutine add_program_requirements_json

    subroutine print_requirement_groups(self, defined_keys)
        class(ui_program), intent(in) :: self
        type(string), optional, intent(in) :: defined_keys(:)
        integer :: i, j, nselected
        if (.not. allocated(self%requirements)) return
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') format_str('INPUT REQUIREMENTS', C_UNDERLINED)
        do i = 1, size(self%requirements)
            nselected = -1
            if (present(defined_keys)) then
                nselected = 0
                do j = 1, size(self%requirements(i)%keys)
                    if (key_is_defined(self%requirements(i)%keys(j), defined_keys)) nselected = nselected + 1
                enddo
            endif
            write(logfhandle,'(a)') trim(self%requirements(i)%label%to_char())//': '// &
                trim(self%requirements(i)%help%to_char())
            write(logfhandle,'(a)') '  Accepted inputs:'
            call print_requirement_inputs(self, self%requirements(i))
            if (nselected >= 0) then
                write(logfhandle,'(a)') '  Supplied: '//int2str(nselected)//'; required: '// &
                    trim(cardinality_text(self%requirements(i)))
            else
                write(logfhandle,'(a)') '  Required: '//trim(cardinality_text(self%requirements(i)))
            endif
        enddo
    end subroutine print_requirement_groups

    logical function has_input_key(self, key)
        class(ui_program), intent(in) :: self
        character(len=*), intent(in) :: key
        has_input_key = list_has_input_key(self%img_ios, key) .or. &
            list_has_input_key(self%file_ios, key) .or. &
            list_has_input_key(self%parm_ios, key) .or. list_has_input_key(self%srch_ctrls, key) .or. &
            list_has_input_key(self%filt_ctrls, key) .or. list_has_input_key(self%mask_ctrls, key) .or. &
            list_has_input_key(self%comp_ctrls, key)
    end function has_input_key

    logical function list_has_input_key(lst, key)
        class(linked_list), intent(in) :: lst
        character(len=*), intent(in) :: key
        type(list_iterator) :: it
        class(*), allocatable :: tmp
        list_has_input_key = .false.
        if (lst%is_empty()) return
        it = lst%begin()
        do while (it%has_value())
            call it%getter(tmp)
            select type(t => tmp)
                type is (ui_program_input)
                    if (t%param%key == trim(key)) then
                        list_has_input_key = .true.
                        if (allocated(tmp)) deallocate(tmp)
                        return
                    endif
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('list_has_input_key: list element is not ui_program_input')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        enddo
    end function list_has_input_key

    logical function key_is_defined(key, defined_keys)
        type(string), intent(in) :: key
        type(string), intent(in) :: defined_keys(:)
        integer :: i
        key_is_defined = .false.
        do i = 1, size(defined_keys)
            if (defined_keys(i) == key) then
                key_is_defined = .true.
                return
            endif
        enddo
    end function key_is_defined

    subroutine print_requirement_inputs(self, requirement)
        class(ui_program), intent(in) :: self
        type(ui_requirement_group), intent(in) :: requirement
        type(chash) :: inputs
        integer :: i
        call inputs%new(size(requirement%keys))
        do i = 1, size(requirement%keys)
            call inputs%push(requirement%keys(i)%to_char(), &
                &input_text_for_key(self, requirement%keys(i)%to_char()))
        enddo
        call inputs%print_key_val_pairs(logfhandle)
        call inputs%kill()
    end subroutine print_requirement_inputs

    function input_text_for_key(self, key) result(text)
        class(ui_program), intent(in) :: self
        character(len=*), intent(in) :: key
        character(len=XLONGSTRLEN) :: text
        logical :: found
        text = ''
        found = .false.
        call find_input_text(self%img_ios, key, text, found)
        if (.not. found) call find_input_text(self%file_ios, key, text, found)
        if (.not. found) call find_input_text(self%parm_ios, key, text, found)
        if (.not. found) call find_input_text(self%srch_ctrls, key, text, found)
        if (.not. found) call find_input_text(self%filt_ctrls, key, text, found)
        if (.not. found) call find_input_text(self%mask_ctrls, key, text, found)
        if (.not. found) call find_input_text(self%comp_ctrls, key, text, found)
        if (.not. found) THROW_HARD('requirement key is not registered: '//trim(key))
    end function input_text_for_key

    subroutine find_input_text(lst, key, text, found)
        class(linked_list), intent(in) :: lst
        character(len=*), intent(in) :: key
        character(len=*), intent(out) :: text
        logical, intent(inout) :: found
        type(list_iterator) :: it
        class(*), allocatable :: item
        if (found .or. lst%is_empty()) return
        it = lst%begin()
        do while (it%has_value())
            call it%getter(item)
            select type(input => item)
                type is(ui_program_input)
                    if (input%param%key%to_char() == trim(key)) then
                        text = input%param%label%to_char()//'; '//input%param%placeholder%to_char()
                        found = .true.
                        if (allocated(item)) deallocate(item)
                        return
                    endif
                class default
                    if (allocated(item)) deallocate(item)
                    THROW_HARD('find_input_text: list element is not ui_program_input')
            end select
            if (allocated(item)) deallocate(item)
            call it%next()
        enddo
    end subroutine find_input_text

    function cardinality_text(requirement) result(text)
        type(ui_requirement_group), intent(in) :: requirement
        character(len=XLONGSTRLEN) :: text
        if (requirement%min_selected == requirement%max_selected) then
            text = 'exactly '//int2str(requirement%min_selected)
        else if (requirement%max_selected == size(requirement%keys)) then
            text = 'at least '//int2str(requirement%min_selected)
        else
            text = int2str(requirement%min_selected)//' to '//int2str(requirement%max_selected)
        endif
    end function cardinality_text

    integer function count_required_in_list( lst ) result( nreq )
        class(linked_list), intent(in) :: lst
        type(list_iterator)            :: it
        class(*), allocatable          :: tmp
        nreq = 0
        if (lst%is_empty()) return
        it = lst%begin()
        do while (it%has_value())
            call it%getter(tmp)
            select type(t => tmp)
                type is (ui_program_input)
                    if (t%param%required) nreq = nreq + 1
                class default
                    THROW_HARD('count_required_in_list: list element is not ui_program_input')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        end do
    end function count_required_in_list

    subroutine append_required_keys_from_list( lst, keys, ireq )
        class(linked_list), intent(in)    :: lst
        type(string),       intent(inout) :: keys(:)
        integer,            intent(inout) :: ireq
        type(list_iterator)   :: it
        class(*), allocatable :: tmp
        if (lst%is_empty()) return
        it = lst%begin()
        do while (it%has_value())
            call it%getter(tmp)
            select type(t => tmp)
            type is (ui_program_input)
                if (t%param%required) then
                    ireq = ireq + 1
                    keys(ireq) = t%param%key
                end if
            class default
                THROW_HARD('append_required_keys_from_list: list element is not ui_program_input')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        end do
    end subroutine append_required_keys_from_list

end module simple_ui_program
