!@descr: module defining the ui_param type, which is used to define input parameters for the simple_ui_program interface
module simple_ui_param
use simple_string, only: string
use simple_error,  only: simple_exception
use simple_ui_descriptor_types, only: ui_choice, ui_choices_are_valid
use simple_string_utils, only: lowercase
implicit none
#include "simple_local_flags.inc"

integer, parameter, public :: UI_PLACEHOLDER_MAX_LEN = 40
public :: ui_placeholder_is_standard

! common input parameter type
type ui_param
    ! deliberately made public (close entaglement with simple_ui_program)
    type(string) :: key
    type(string) :: keytype
    type(string) :: label
    type(string) :: help
    type(string) :: placeholder
    type(string) :: units
    type(string) :: example_value
    type(string) :: cval_default
    real         :: rval_default = 0.
    logical      :: required = .true.
    logical      :: has_default = .false.
    type(ui_choice), allocatable :: choices(:)
contains
    procedure, private :: set_param_1
    procedure, private :: set_param_2
    generic            :: set_param => set_param_1, set_param_2
    procedure          :: set_choices
    procedure          :: set_placeholder
    procedure          :: set_generated_default
    final              :: finalize
end type ui_param

contains

    subroutine set_param_1( self, key, keytype, label, help, placeholder, required, default_value, choices )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        real,                   intent(in)    :: default_value
        type(ui_choice), optional, intent(in)  :: choices(:)
        self%key         = trim(key)
        self%keytype     = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%required    = required
        self%has_default = .not. self%required
        if( self%has_default ) self%rval_default = default_value
        call self%set_choices(choices)
        call self%set_placeholder(placeholder, '')
    end subroutine set_param_1

    subroutine set_param_2( self, key, keytype, label, help, placeholder, required, default_value, choices )
        class(ui_param), intent(inout) :: self
        character(len=*),       intent(in)    :: key, keytype, label, help, placeholder
        logical,                intent(in)    :: required
        character(len=*),       intent(in)    :: default_value
        type(ui_choice), optional, intent(in)  :: choices(:)
        self%key         = trim(key)
        self%keytype     = trim(keytype)
        self%label       = trim(label)
        self%help        = trim(help)
        self%required    = required
        self%has_default = .not. self%required
        if( self%has_default ) self%cval_default = trim(default_value)
        call self%set_choices(choices)
        call self%set_placeholder(placeholder, '')
        if( self%has_default .and. (self%keytype%to_char() == 'binary' .or. &
            &self%keytype%to_char() == 'multi') )then
            if( .not. choice_is_available(self%choices, self%cval_default%to_char()) )then
                THROW_HARD('UI parameter default is not a declared choice: '//self%key%to_char())
            endif
        endif
    end subroutine set_param_2

    subroutine set_choices( self, choices )
        class(ui_param), intent(inout) :: self
        type(ui_choice), optional, intent(in) :: choices(:)
        if( allocated(self%choices) ) deallocate(self%choices)
        if( present(choices) )then
            allocate(self%choices(size(choices)))
            self%choices = choices
        else
            allocate(self%choices(0))
        endif
        if( .not. ui_choices_are_valid(self%keytype%to_char(), self%choices) )then
            THROW_HARD('Invalid choices for UI parameter: '//self%key%to_char())
        endif
    end subroutine set_choices

    subroutine set_placeholder( self, placeholder, example_value )
        class(ui_param), intent(inout) :: self
        character(len=*), intent(in)   :: placeholder, example_value
        character(len=:), allocatable  :: keytype
        character(len=UI_PLACEHOLDER_MAX_LEN) :: standard_placeholder
        keytype = trim(self%keytype%to_char())
        self%units = infer_units(keytype, self%label%to_char(), self%help%to_char(), placeholder)
        self%example_value = trim(example_value)
        standard_placeholder = standardize_placeholder(self%key%to_char(), keytype, self%units%to_char(), &
            &self%example_value%to_char())
        select case( keytype )
            case('file')
                if( trim(standard_placeholder) /= 'e.g. input.file' )then
                    self%placeholder = trim(standard_placeholder)
                else if( len_trim(placeholder) > 0 )then
                    self%placeholder = trim(placeholder)
                else
                    self%placeholder = trim(standard_placeholder)
                endif
            case('dir')
                if( len_trim(placeholder) > 0 )then
                    self%placeholder = trim(placeholder)
                else
                    self%placeholder = trim(standard_placeholder)
                endif
            case default
                self%placeholder = trim(standard_placeholder)
        end select
        if( len_trim(self%placeholder%to_char()) > UI_PLACEHOLDER_MAX_LEN )then
            THROW_HARD('UI parameter placeholder exceeds 40 characters: '//self%key%to_char())
        endif
    end subroutine set_placeholder

    subroutine set_generated_default( self, value, found )
        class(ui_param), intent(inout) :: self
        character(len=*), intent(in)   :: value
        logical,          intent(in)   :: found
        integer :: ios
        logical :: had_default
        if( .not. found ) return
        had_default = self%has_default
        select case( trim(self%keytype%to_char()) )
            case('num', 'int', 'float')
                read(value, *, iostat=ios) self%rval_default
                if( ios /= 0 )then
                    THROW_HARD('Generated UI default is not numeric: '//self%key%to_char())
                endif
            case default
                if( self%keytype%to_char() == 'binary' .or. self%keytype%to_char() == 'multi' )then
                    if( .not. choice_is_available(self%choices, trim(value)) )then
                        ! A baseline parameter value may be valid globally but outside
                        ! this program's declared subset. Retain its validated local
                        ! default; an input without one remains a descriptor error.
                        if( had_default ) return
                        THROW_HARD('Generated UI default is not a declared choice: '//self%key%to_char())
                    endif
                endif
                self%cval_default = trim(value)
        end select
        self%has_default = .true.
        if( self%keytype%to_char() == 'num' .or. self%keytype%to_char() == 'int' .or. &
            &self%keytype%to_char() == 'float' )then
            call self%set_placeholder(self%placeholder%to_char(), value)
        endif
    end subroutine set_generated_default

    pure logical function choice_is_available( choices, value )
        type(ui_choice),  intent(in) :: choices(:)
        character(len=*), intent(in) :: value
        integer :: i
        choice_is_available = .false.
        do i = 1, size(choices)
            if( choices(i)%value%to_char() == trim(value) )then
                choice_is_available = .true.
                return
            endif
        enddo
    end function choice_is_available

    pure function standardize_placeholder( key, keytype, units, example_value ) result( placeholder )
        character(len=*), intent(in) :: key, keytype, units, example_value
        character(len=UI_PLACEHOLDER_MAX_LEN) :: placeholder
        ! Project files are a semantic subset of generic files.  Keep their
        ! format visible in every UI context rather than suggesting an MRC.
        if( index(trim(key), 'projfile') == 1 )then
            placeholder = 'e.g. input.simple'
            return
        endif
        if( trim(key) == 'subprojname' )then
            placeholder = 'e.g. myproject.simple'
            return
        endif
        if( trim(key) == 'projtab' )then
            placeholder = 'e.g. projtab.txt'
            return
        endif
        select case( trim(key) )
            case('blocktree')
                placeholder = 'e.g. pool_block_tree.bin'
                return
            case('boxfile')
                placeholder = 'e.g. coords.box'
                return
            case('boxtab')
                placeholder = 'e.g. boxes.txt'
                return
            case('ciffile')
                placeholder = 'e.g. molecule.cif'
                return
            case('deftab')
                placeholder = 'e.g. deftab.txt'
                return
            case('filetab')
                placeholder = 'e.g. input_files.txt'
                return
            case('fsc')
                placeholder = 'e.g. fsc.bin'
                return
            case('gainref')
                placeholder = 'e.g. gainref.mrc'
                return
            case('oritab')
                placeholder = 'e.g. oritab.txt'
                return
            case('oritab2')
                placeholder = 'e.g. oritab2.txt'
                return
            case('outstk')
                placeholder = 'e.g. output.mrcs'
                return
            case('outvol')
                placeholder = 'e.g. volume.mrc'
                return
            case('pdbfile', 'pdbfile2')
                placeholder = 'e.g. molecule.pdb'
                return
            case('pdbfiles')
                placeholder = 'e.g. pdb_files.txt'
                return
            case('pdbout')
                placeholder = 'e.g. output.pdb'
                return
            case('pickrefs')
                placeholder = 'e.g. pickrefs.mrc'
                return
            case('plaintexttab')
                placeholder = 'e.g. params.txt'
                return
            case('refs')
                placeholder = 'e.g. references.mrc'
                return
            case('rmsd_file')
                placeholder = 'e.g. atom_rmsd.bin'
                return
            case('star_mic')
                placeholder = 'e.g. micrographs.star'
                return
            case('star_model')
                placeholder = 'e.g. model.star'
                return
            case('star_ptcl')
                placeholder = 'e.g. particles.star'
                return
            case('starfile')
                placeholder = 'e.g. input.star'
                return
            case('stk')
                placeholder = 'e.g. particles.mrcs'
                return
            case('stk2')
                placeholder = 'e.g. particles2.mrcs'
                return
            case('stk_backgr')
                placeholder = 'e.g. background_pspec.mrcs'
                return
            case('stk_den')
                placeholder = 'e.g. denoised.mrcs'
                return
            case('stk_traj')
                placeholder = 'e.g. trajectory.mrcs'
                return
            case('stktab')
                placeholder = 'e.g. stktab.txt'
                return
            case('stktab_den')
                placeholder = 'e.g. stktab_den.txt'
                return
            case('vol1')
                placeholder = 'e.g. volume.mrc'
                return
            case('vol2')
                placeholder = 'e.g. volume2.mrc'
                return
            case('vol3')
                placeholder = 'e.g. volume3.mrc'
                return
            case('vol_odd')
                placeholder = 'e.g. odd.mrc'
                return
            case('vol_even')
                placeholder = 'e.g. even.mrc'
                return
        end select
        select case( trim(keytype) )
            case('binary', 'multi', 'hidden_dir', 'hidden_float', 'hidden_int')
                placeholder = ''
            case('dir')
                placeholder = 'e.g. /path/to/folder'
            case('file')
                placeholder = 'e.g. input.file'
            case('num', 'int', 'float')
                if( len_trim(example_value) > 0 )then
                    placeholder = 'e.g. '//trim(normalize_numeric_example(example_value))
                else
                    select case( trim(key) )
                        case('mskdiam', 'mskdiam_detect')
                            placeholder = 'e.g. 180'
                        case default
                            placeholder = 'e.g. 10'
                    end select
                endif
                if( len_trim(units) > 0 ) placeholder = trim(placeholder)//' '//trim(units)
            case('str', 'string')
                placeholder = 'e.g. value'
            case default
                placeholder = 'e.g. value'
        end select
    end function standardize_placeholder

    pure logical function ui_placeholder_is_standard( key, keytype, placeholder, units, example_value )
        character(len=*), intent(in) :: key, keytype, placeholder, units, example_value
        select case( trim(keytype) )
            case('file', 'dir')
                ui_placeholder_is_standard = len_trim(placeholder) > 0 .and. &
                    &len_trim(placeholder) <= UI_PLACEHOLDER_MAX_LEN
            case default
                ui_placeholder_is_standard = trim(placeholder) == &
                    &trim(standardize_placeholder(key, keytype, units, example_value))
        end select
    end function ui_placeholder_is_standard

    pure function normalize_numeric_example( value ) result( normalized )
        character(len=*), intent(in) :: value
        character(len=40)            :: normalized
        normalized = trim(value)
        if( index(normalized, 'e') == 0 .and. index(normalized, 'E') == 0 )then
            if( normalized(len_trim(normalized):len_trim(normalized)) == '.' )then
                normalized = normalized(:len_trim(normalized) - 1)
            endif
        endif
    end function normalize_numeric_example

    pure function infer_units( keytype, label, help, placeholder ) result( units )
        character(len=*), intent(in) :: keytype, label, help, placeholder
        character(len=16)            :: units
        character(len=:), allocatable :: text
        units = ''
        select case( trim(keytype) )
            case('num', 'int', 'float')
            case default
                return
        end select
        text = lowercase(trim(label)//' '//trim(help)//' '//trim(placeholder))
        if( index(text, 'e/a2') > 0 .or. index(text, 'e/angstrom') > 0 )then
            units = 'e/Angstrom^2'
        else if( index(text, 'angstrom^2') > 0 .or. index(text, 'angstroms^2') > 0 )then
            units = 'Angstroms^2'
        else if( index(text, 'angstrom') > 0 .or. index(text, ' angs') > 0 .or. &
            &index(text, '(in a)') > 0 .or. index(text, ' in a {') > 0 )then
            units = 'Angstroms'
        else if( index(text, 'pixel') > 0 )then
            units = 'pixels'
        else if( index(text, 'degree') > 0 )then
            units = 'degrees'
        else if( index(text, 'radian') > 0 )then
            units = 'radians'
        else if( index(text, 'micron') > 0 )then
            units = 'microns'
        else if( index(text, 'kilovolt') > 0 .or. index(text, '(kv)') > 0 .or. index(text, ' in kv') > 0 )then
            units = 'kV'
        else if( index(text, 'millimet') > 0 .or. index(text, '(mm)') > 0 .or. index(text, 'in mm') > 0 )then
            units = 'mm'
        else if( index(text, 'second') > 0 )then
            units = 'seconds'
        else if( index(text, 'minute') > 0 )then
            units = 'minutes'
        else if( index(text, 'hour') > 0 )then
            units = 'hours'
        endif
    end function infer_units

    subroutine finalize(self)
        type(ui_param), intent(inout) :: self
        call self%key%kill()
        call self%keytype%kill()
        call self%label%kill()
        call self%help%kill()
        call self%placeholder%kill()
        call self%units%kill()
        call self%cval_default%kill()
    end subroutine finalize

end module simple_ui_param
