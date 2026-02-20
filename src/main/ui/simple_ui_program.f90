!@descr: module defining the main user interface program type and associated methods
module simple_ui_program
use simple_core_module_api
use simple_ansi_ctrls
use simple_linked_list, only: linked_list, list_iterator
use simple_ui_param,    only: ui_param
implicit none
#include "simple_local_flags.inc"

integer, parameter :: UI_IMG=1, UI_PARM=2, UI_ALT=3, UI_SRCH=4, UI_FILT=5, UI_MASK=6, UI_COMP=7

! production-level program interface for simple_exec, single_exec & simple_stream executables
type :: ui_program
    type(string) :: name
    type(string) :: descr_short
    type(string) :: descr_long
    type(string) :: executable
    type(string) :: gui_submenu_list
    logical      :: advanced = .true.
    ! image input/output
    type(linked_list) :: img_ios
    ! parameter input/output
    type(linked_list) :: parm_ios
    ! alternative inputs
    type(linked_list) :: alt_ios
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
    procedure          :: print_prg_descr_long
    procedure          :: write2json
    procedure          :: get_name
    procedure          :: get_executable
    procedure          :: get_nrequired_keys
    procedure          :: get_required_keys
    procedure          :: requires_sp_project
    procedure, private :: kill
end type ui_program

contains

    subroutine new( self, name, descr_short, descr_long, executable, sp_required, gui_advanced, gui_submenu_list )
        class(ui_program),          intent(inout) :: self
        character(len=*),           intent(in)    :: name, descr_short, descr_long, executable
        logical,                    intent(in)    :: sp_required
        logical,          optional, intent(in)    :: gui_advanced
        character(len=*), optional, intent(in)    :: gui_submenu_list
        call self%kill
        self%name        = trim(name)
        self%descr_short = trim(descr_short)
        self%descr_long  = trim(descr_long)
        self%executable  = trim(executable)
        self%sp_required = sp_required
        self%exists      = .true.
        if(present(gui_advanced)    ) self%advanced         = gui_advanced
        if(present(gui_submenu_list)) self%gui_submenu_list = gui_submenu_list
    end subroutine new

    subroutine add_input_num( self, which, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value, &
                            gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        character(len=*),           intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                    intent(in)    :: required
        real,                       intent(in)    :: default_value
        character(len=*), optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,          optional, intent(in)    :: gui_advanced, gui_online
        type(ui_param) :: p
        call p%set_param(key, keytype, descr_short, descr_long, descr_placeholder, required, default_value)
        call p%apply_gui_overrides(gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online)
        select case (which)
            case (UI_IMG);  call self%img_ios%push_back(p)
            case (UI_PARM); call self%parm_ios%push_back(p)
            case (UI_ALT);  call self%alt_ios%push_back(p)
            case (UI_SRCH); call self%srch_ctrls%push_back(p)
            case (UI_FILT); call self%filt_ctrls%push_back(p)
            case (UI_MASK); call self%mask_ctrls%push_back(p)
            case (UI_COMP); call self%comp_ctrls%push_back(p)
            case default
                THROW_HARD('which field selector: '//int2str(which)//' is unsupported; ui_program::add_input_num')
        end select
    end subroutine add_input_num

    subroutine add_input_str( self, which, key, keytype, descr_short, descr_long, descr_placeholder, required, default_value, &
                            gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        character(len=*),           intent(in)    :: key, keytype, descr_short, descr_long, descr_placeholder
        logical,                    intent(in)    :: required
        character(len=*),           intent(in)    :: default_value
        character(len=*), optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,          optional, intent(in)    :: gui_advanced, gui_online
        type(ui_param) :: p
        call p%set_param(key, keytype, descr_short, descr_long, descr_placeholder, required, default_value)
        call p%apply_gui_overrides(gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online)
        select case (which)
            case (UI_IMG);  call self%img_ios%push_back(p)
            case (UI_PARM); call self%parm_ios%push_back(p)
            case (UI_ALT);  call self%alt_ios%push_back(p)
            case (UI_SRCH); call self%srch_ctrls%push_back(p)
            case (UI_FILT); call self%filt_ctrls%push_back(p)
            case (UI_MASK); call self%mask_ctrls%push_back(p)
            case (UI_COMP); call self%comp_ctrls%push_back(p)
            case default
                THROW_HARD('which field selector: '//int2str(which)//' is unsupported; ui_program::add_input_str')
        end select
    end subroutine add_input_str

    subroutine add_input_param( self, which, param, descr_short_override, descr_long_override, descr_placeholder_override,&
    &required_override, gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online )
        class(ui_program),          intent(inout) :: self
        integer,                    intent(in)    :: which
        type(ui_param),             intent(in)    :: param
        character(len=*), optional, intent(in)    :: descr_short_override, descr_long_override, descr_placeholder_override
        logical,          optional, intent(in)    :: required_override
        character(len=*), optional, intent(in)    :: gui_submenu, gui_exclusive_group, gui_active_flags
        logical,          optional, intent(in)    :: gui_advanced, gui_online
        type(ui_param) :: p
        p = param
        if( present(descr_short_override)       ) p%descr_short       = descr_short_override
        if( present(descr_long_override)        ) p%descr_long        = descr_long_override
        if( present(descr_placeholder_override) ) p%descr_placeholder = descr_placeholder_override
        if( present(required_override)          ) p%required          = required_override
        call p%apply_gui_overrides(gui_submenu, gui_exclusive_group, gui_active_flags, gui_advanced, gui_online)
        select case (which)
            case (UI_IMG);  call self%img_ios%push_back(p)
            case (UI_PARM); call self%parm_ios%push_back(p)
            case (UI_ALT);  call self%alt_ios%push_back(p)
            case (UI_SRCH); call self%srch_ctrls%push_back(p)
            case (UI_FILT); call self%filt_ctrls%push_back(p)
            case (UI_MASK); call self%mask_ctrls%push_back(p)
            case (UI_COMP); call self%comp_ctrls%push_back(p)
            case default
                THROW_HARD('which field selector: '//int2str(which)//' is unsupported; ui_program::add_input_param')
        end select
    end subroutine add_input_param

    subroutine print_ui( self )
        class(ui_program), intent(in) :: self
        type(chash) :: ch
        write(logfhandle,'(a)') ''
        write(logfhandle,'(a)') '>>> PROGRAM INFO'
        call ch%new(4)
        call ch%push('name',        self%name%to_char())
        call ch%push('descr_short', self%descr_short%to_char())
        call ch%push('descr_long',  self%descr_long%to_char())
        call ch%push('executable',  self%executable%to_char())
        call ch%print_key_val_pairs(logfhandle)
        call ch%kill
        if (.not. self%img_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT',     C_UNDERLINED)
            call print_param_list(self%img_ios)
        endif
        if (.not. self%parm_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
            call print_param_list(self%parm_ios)
        endif
        if (.not. self%alt_ios%is_empty()) then
            write(logfhandle,'(a)') ''
            write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS',     C_UNDERLINED)
            call print_param_list(self%alt_ios)
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

    subroutine print_cmdline( self )
        class(ui_program), intent(in) :: self
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
        write(logfhandle,'(a)') 'Required input parameters in ' // format_str('bold', C_BOLD) // ' (ensure terminal support)'
        if (.not. self%img_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('IMAGE INPUT/OUTPUT', C_UNDERLINED)
            call print_param_hash(self%img_ios)
        end if
        if (.not. self%parm_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('PARAMETER INPUT/OUTPUT', C_UNDERLINED)
            call print_param_hash(self%parm_ios)
        end if
        if (.not. self%alt_ios%is_empty()) then
            write(logfhandle,'(a)') format_str('ALTERNATIVE INPUTS', C_UNDERLINED)
            call print_param_hash(self%alt_ios)
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
    end subroutine print_cmdline

    subroutine print_prg_descr_long( self )
        class(ui_program), intent(in) :: self
        write(logfhandle,'(a)') self%descr_long%to_char()
    end subroutine print_prg_descr_long

    subroutine write2json( self )
        use json_module
        class(ui_program), intent(in) :: self
        type(json_core)               :: json
        type(json_value), pointer     :: program_entry, program
        ! JSON init
        call json%initialize()
        call json%create_object(program_entry,'')
        call json%create_object(program, self%name%to_char())
        call json%add(program_entry, program)
        ! program section
        call json%add(program, 'name',        self%name%to_char())
        call json%add(program, 'descr_short', self%descr_short%to_char())
        call json%add(program, 'descr_long',  self%descr_long%to_char())
        call json%add(program, 'executable',  self%executable%to_char())
        ! all sections (linked lists)
        call create_section_from_list(json, program_entry, 'image input/output',     self%img_ios)
        call create_section_from_list(json, program_entry, 'parameter input/output', self%parm_ios)
        call create_section_from_list(json, program_entry, 'alternative inputs',     self%alt_ios)
        call create_section_from_list(json, program_entry, 'search controls',        self%srch_ctrls)
        call create_section_from_list(json, program_entry, 'filter controls',        self%filt_ctrls)
        call create_section_from_list(json, program_entry, 'mask controls',          self%mask_ctrls)
        call create_section_from_list(json, program_entry, 'computer controls',      self%comp_ctrls)
        ! write & clean
        call json%print(program_entry, self%name%to_char()//'.json')
        if( json%failed() )then
            THROW_HARD('json input/output error for program: '//self%name%to_char())
        endif
        call json%destroy(program_entry)
    end subroutine write2json

    function get_name( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%name
    end function get_name

    function get_executable( self ) result( name )
        class(ui_program), intent(in) :: self
        type(string) :: name
        name = self%executable
    end function get_executable

    integer function get_nrequired_keys( self )
        class(ui_program), intent(in) :: self
        get_nrequired_keys = count_required_in_list(self%img_ios)   + &
                            count_required_in_list(self%parm_ios)   + &
                            count_required_in_list(self%srch_ctrls) + &
                            count_required_in_list(self%filt_ctrls) + &
                            count_required_in_list(self%mask_ctrls) + &
                            count_required_in_list(self%comp_ctrls)
        ! Preserve legacy behavior: if *no required keys* exist anywhere above,
        ! but there are alternative inputs, then require at least one alt input.
        if (get_nrequired_keys == 0 .and. .not. self%alt_ios%is_empty()) then
            get_nrequired_keys = 1
        end if
    end function get_nrequired_keys

    function get_required_keys( self ) result( keys )
        class(ui_program), intent(in) :: self
        type(string), allocatable     :: keys(:)
        integer                       :: nreq, ireq
        type(list_iterator)           :: it
        class(*), allocatable         :: tmp
        ! count # required
        nreq = self%get_nrequired_keys()
        if (nreq <= 0) return
        allocate(keys(nreq))
        ireq = 0
        call append_required_keys_from_list(self%img_ios,    keys, ireq)
        call append_required_keys_from_list(self%parm_ios,   keys, ireq)
        call append_required_keys_from_list(self%alt_ios,    keys, ireq)
        call append_required_keys_from_list(self%srch_ctrls, keys, ireq)
        call append_required_keys_from_list(self%filt_ctrls, keys, ireq)
        call append_required_keys_from_list(self%mask_ctrls, keys, ireq)
        call append_required_keys_from_list(self%comp_ctrls, keys, ireq)
        ! legacy alt_ios fallback: require at least one alt key if nothing else required
        if (ireq == 0 .and. nreq == 1 .and. .not. self%alt_ios%is_empty()) then
            it = self%alt_ios%begin()
            call it%getter(tmp)
            select type(t => tmp)
                type is (ui_param)
                    ireq = 1
                    keys(1) = t%key
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('get_required_keys: alt_ios element is not ui_param')
            end select
            if (allocated(tmp)) deallocate(tmp)
        end if
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
        call self%descr_short%kill()
        call self%descr_long%kill()
        call self%executable%kill()
        call self%gui_submenu_list%kill()
        call self%img_ios%kill()
        call self%parm_ios%kill()
        call self%alt_ios%kill()
        call self%srch_ctrls%kill()
        call self%filt_ctrls%kill()
        call self%mask_ctrls%kill()
        call self%comp_ctrls%kill()
        self%advanced    = .true.
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
                type is (ui_param)
                    write(logfhandle,'(a,1x,i3)') '>>> PARAMETER #', i
                    call ch%new(6)
                    call ch%push('key',               t%key%to_char())
                    call ch%push('keytype',           t%keytype%to_char())
                    call ch%push('descr_short',       t%descr_short%to_char())
                    call ch%push('descr_long',        t%descr_long%to_char())
                    call ch%push('descr_placeholder', t%descr_placeholder%to_char())
                    call ch%push('required', merge('T','F', t%required))
                    call ch%print_key_val_pairs(logfhandle)
                    call ch%kill
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('print_param_list: list element is not ui_param')
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
            type is (ui_param)
                call ch%push(t%key%to_char(), t%descr_short%to_char()//'; '//t%descr_placeholder%to_char())
                keys(i) = t%key%to_char()
                req(i)  = t%required
            class default
                if (allocated(tmp)) deallocate(tmp)
                THROW_HARD('print_param_hash: list element is not ui_param')
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
        character(len=STDLEN)     :: options_str, before
        character(len=KEYLEN)     :: args(8)
        integer                   :: j, nargs
        logical                   :: found, param_is_multi, param_is_binary, exception
        call json%create_array(section, trim(name))
        if (.not. lst%is_empty()) then
            it = lst%begin()
            do while (it%has_value())
                call it%getter(tmp)
                select type(t => tmp)
                type is (ui_param)
                    call json%create_object(entry, t%key%to_char())
                    call json%add(entry, 'key',               t%key%to_char())
                    call json%add(entry, 'keytype',           t%keytype%to_char())
                    call json%add(entry, 'descr_short',       t%descr_short%to_char())
                    call json%add(entry, 'descr_long',        t%descr_long%to_char())
                    call json%add(entry, 'descr_placeholder', t%descr_placeholder%to_char())
                    call json%add(entry, 'required',          t%required)
                    ! Optional: emit defaults when not required
                    if (.not. t%required) then
                        if (t%keytype%to_char() == 'num') then
                            call json%add(entry, 'default', real(t%rval_default,dp))
                        else
                            call json%add(entry, 'default', t%cval_default%to_char())
                        end if
                    end if
                    ! Optional: emit GUI fields (only if you want them in JSON)
                    if (len_trim(t%gui_submenu%to_char()) > 0) then
                        call json%add(entry, 'gui_submenu', t%gui_submenu%to_char())
                    end if
                    if (len_trim(t%exclusive_group%to_char()) > 0) then
                        call json%add(entry, 'exclusive_group', t%exclusive_group%to_char())
                    end if
                    if (len_trim(t%active_flags%to_char()) > 0) then
                        call json%add(entry, 'active_flags', t%active_flags%to_char())
                    end if
                    call json%add(entry, 'advanced', t%advanced)
                    call json%add(entry, 'online',   t%online)
                    ! options parsing (multi/binary)
                    param_is_multi  = (t%keytype%to_char() == 'multi')
                    param_is_binary = (t%keytype%to_char() == 'binary')
                    if (param_is_multi .or. param_is_binary) then
                        options_str = t%descr_placeholder%to_char()
                        call split(options_str, '(', before)
                        call split(options_str, ')', before)
                        call parsestr(before, '|', args, nargs)
                        exception = (param_is_binary .and. nargs /= 2) .or. (param_is_multi .and. nargs < 2)
                        if (exception) then
                            write(logfhandle,*) 'Poorly formatted options string for entry ', t%key%to_char()
                            write(logfhandle,*) t%descr_placeholder%to_char()
                            THROW_HARD('Bad options string formatting in UI JSON export')
                        end if
                        call json%add(entry, 'options', args(1:nargs))
                        do j = 1, nargs
                            call json%update(entry, 'options['//int2str(j)//']', trim(args(j)), found)
                        end do
                    end if
                    call json%add(section, entry)
                class default
                    if (allocated(tmp)) deallocate(tmp)
                    THROW_HARD('create_section_from_list: list element is not ui_param')
                end select
                if (allocated(tmp)) deallocate(tmp)
                call it%next()
            end do
        end if
        call json%add(program_entry, section)
    end subroutine create_section_from_list

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
                type is (ui_param)
                    if (t%required) nreq = nreq + 1
                class default
                    THROW_HARD('count_required_in_list: list element is not ui_param')
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
            type is (ui_param)
                if (t%required) then
                    ireq = ireq + 1
                    keys(ireq) = t%key
                end if
            class default
                THROW_HARD('append_required_keys_from_list: list element is not ui_param')
            end select
            if (allocated(tmp)) deallocate(tmp)
            call it%next()
        end do
    end subroutine append_required_keys_from_list

end module simple_ui_program
