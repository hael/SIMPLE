! View file: top-level organization for a refactored simple_parameters module.
! This is intentionally architectural rather than compilable.

module simple_parameters_refactor_view
use simple_core_module_api
use simple_cmdline, only: cmdline
use simple_ui_program, only: ui_program
use simple_ui, only: get_prg_ptr
use simple_parameters_proto_registry, only: param_registry
implicit none

public :: parameters
private

type :: parameters
    ! The full flat field inventory remains in the same overall style as today.
    ! See simple_parameters_type_reference.f90 for the real declaration surface.
    !
    ! The key design change is not nested storage. It is:
    ! 1. defaults live on the declarations
    ! 2. input parsing is registry-driven
    ! 3. semantic work is split into named phases
    type(ui_program), pointer :: ptr2prg => null()

    ! [same flat field inventory as simple_parameters_type_reference.f90]

contains
    procedure :: new
    procedure :: is_final_planned_iter
    procedure, private :: set_img_format

    ! Reset and parsing
    procedure, private :: reset_defaults
    procedure, private :: parse_inputs
    procedure, private :: bind_input_char3
    procedure, private :: bind_input_character
    procedure, private :: bind_input_string
    procedure, private :: bind_input_integer
    procedure, private :: bind_input_real

    ! Execution/bootstrap
    procedure, private :: init_execution_context
    procedure, private :: attach_program_context
    procedure, private :: prepare_execution_directory

    ! Project/input resolution
    procedure, private :: normalize_path_inputs
    procedure, private :: resolve_project_requirements
    procedure, private :: resolve_project_inputs
    procedure, private :: resolve_stack_and_reference_inputs
    procedure, private :: resolve_volume_inputs

    ! Parallel/runtime policy
    procedure, private :: derive_partitioning_defaults
    procedure, private :: derive_fractional_update_flags
    procedure, private :: derive_threading_settings
    procedure, private :: validate_walltime

    ! Image/Fourier/scaling
    procedure, private :: derive_movie_scaling
    procedure, private :: derive_box_and_sampling
    procedure, private :: derive_mask_defaults
    procedure, private :: derive_frequency_defaults
    procedure, private :: derive_filter_mode_flags
    procedure, private :: derive_refine_mode_flags
    procedure, private :: derive_misc_image_flags

    ! Validation
    procedure, private :: validate_enum_values
    procedure, private :: validate_file_combinations
    procedure, private :: validate_dimensions
    procedure, private :: validate_image_settings
    procedure, private :: validate_refine_settings
end type parameters

interface
    module subroutine new(self, cline, silent)
        class(parameters), intent(inout) :: self
        class(cmdline),    intent(inout) :: cline
        logical, optional, intent(in)    :: silent
    end subroutine new
end interface

contains

    logical function is_final_planned_iter(self)
        class(parameters), intent(in) :: self
        is_final_planned_iter = .false.
    end function is_final_planned_iter

    ! Constructor flow:
    !
    ! module subroutine new(self, cline, silent)
    !     self = parameters()
    !     call self%parse_inputs(cline)
    !
    !     call self%init_execution_context(cline)
    !     call self%attach_program_context(cline)
    !     call self%prepare_execution_directory(cline)
    !
    !     call self%normalize_path_inputs(cline)
    !     call self%resolve_project_requirements(cline)
    !     call self%resolve_project_inputs(cline)
    !     call self%resolve_stack_and_reference_inputs(cline)
    !     call self%resolve_volume_inputs(cline)
    !
    !     call self%derive_partitioning_defaults(cline)
    !     call self%derive_fractional_update_flags(cline)
    !     call self%derive_threading_settings(cline)
    !     call self%validate_walltime(cline)
    !
    !     call self%derive_movie_scaling(cline)
    !     call self%derive_box_and_sampling(cline)
    !     call self%derive_mask_defaults(cline)
    !     call self%derive_frequency_defaults(cline)
    !     call self%derive_filter_mode_flags(cline)
    !     call self%derive_refine_mode_flags(cline)
    !     call self%derive_misc_image_flags(cline)
    !
    !     call self%validate_enum_values(cline)
    !     call self%validate_file_combinations(cline)
    !     call self%validate_dimensions(cline)
    !     call self%validate_image_settings(cline)
    !     call self%validate_refine_settings(cline)
    ! end subroutine new

end module simple_parameters_refactor_view
