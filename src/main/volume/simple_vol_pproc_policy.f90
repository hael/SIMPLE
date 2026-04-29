!@descr: per-state automask and nonuniform-filter policy decisions consumed by Cartesian assembly
module simple_vol_pproc_policy
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

private

public :: vol_pproc_plan
public :: plan_state_postprocess
public :: state_mask_is_compatible
public :: AUTOMASK_ACTION_NONE
public :: AUTOMASK_ACTION_REGENERATE
public :: AUTOMASK_ACTION_REUSE
public :: NU_MASK_SOURCE_NONE
public :: NU_MASK_SOURCE_FRESH_AUTOMASK
public :: NU_MASK_SOURCE_EXISTING_AUTOMASK
public :: NU_MASK_SOURCE_SPHERICAL

integer, parameter :: AUTOMASK_ACTION_NONE       = 0
integer, parameter :: AUTOMASK_ACTION_REGENERATE = 1
integer, parameter :: AUTOMASK_ACTION_REUSE       = 2

integer, parameter :: NU_MASK_SOURCE_NONE              = 0
integer, parameter :: NU_MASK_SOURCE_FRESH_AUTOMASK    = 1
integer, parameter :: NU_MASK_SOURCE_EXISTING_AUTOMASK = 2
integer, parameter :: NU_MASK_SOURCE_SPHERICAL         = 3

type :: vol_pproc_plan
    type(string) :: mskfile_state
    integer      :: automask_action              = AUTOMASK_ACTION_NONE
    integer      :: nu_mask_source               = NU_MASK_SOURCE_NONE
    logical      :: l_automask_enabled            = .false.
    logical      :: automask_tight                = .false.
    logical      :: l_state_mask_exists           = .false.
    logical      :: l_state_mask_compatible       = .false.
    logical      :: l_state_mask_incompatible     = .false.
end type vol_pproc_plan

contains

    subroutine plan_state_postprocess( params, state, which_iter, l_nonuniform_mode, plan )
        class(parameters),            intent(in)    :: params
        integer,                      intent(in)    :: state
        integer,                      intent(in)    :: which_iter
        logical,                      intent(in)    :: l_nonuniform_mode
        type(vol_pproc_plan), intent(inout) :: plan
        plan%mskfile_state              = string(AUTOMASK_FBODY)//int2str_pad(state,2)//string(MRC_EXT)
        plan%automask_action            = AUTOMASK_ACTION_NONE
        plan%nu_mask_source             = NU_MASK_SOURCE_NONE
        plan%l_automask_enabled         = trim(params%automsk).ne.'no'
        plan%automask_tight             = trim(params%automsk).eq.'tight'
        plan%l_state_mask_exists        = .false.
        plan%l_state_mask_compatible    = .false.
        plan%l_state_mask_incompatible  = .false.
        if( .not. plan%l_automask_enabled )then
            if( l_nonuniform_mode ) plan%nu_mask_source = NU_MASK_SOURCE_SPHERICAL
            return
        endif
        call state_mask_is_compatible(plan%mskfile_state, params%box_crop, params%smpd_crop, &
            &plan%l_state_mask_exists, plan%l_state_mask_compatible)
        plan%l_state_mask_incompatible = plan%l_state_mask_exists .and. (.not. plan%l_state_mask_compatible)
        if( should_regenerate_automask(params, which_iter, plan%l_state_mask_exists, &
            &plan%l_state_mask_compatible) )then
            plan%automask_action = AUTOMASK_ACTION_REGENERATE
        else if( plan%l_state_mask_compatible )then
            plan%automask_action = AUTOMASK_ACTION_REUSE
        endif
        if( l_nonuniform_mode )then
            select case( plan%automask_action )
                case( AUTOMASK_ACTION_REGENERATE )
                    plan%nu_mask_source = NU_MASK_SOURCE_FRESH_AUTOMASK
                case( AUTOMASK_ACTION_REUSE )
                    plan%nu_mask_source = NU_MASK_SOURCE_EXISTING_AUTOMASK
                case default
                    plan%nu_mask_source = NU_MASK_SOURCE_SPHERICAL
            end select
        endif
    end subroutine plan_state_postprocess

    subroutine state_mask_is_compatible( mskfile_state, box, smpd, exists, compatible )
        class(string), intent(in)  :: mskfile_state
        integer,       intent(in)  :: box
        real,          intent(in)  :: smpd
        logical,       intent(out) :: exists
        logical,       intent(out) :: compatible
        real    :: smpd_mask
        integer :: ldim_mask(3), nptcls_mask
        exists      = file_exists(mskfile_state)
        compatible  = .false.
        if( .not. exists ) return
        call find_ldim_nptcls(mskfile_state, ldim_mask, nptcls_mask)
        smpd_mask = find_img_smpd(mskfile_state)
        compatible = ldim_mask(1) == box .and. ldim_mask(2) == box .and. ldim_mask(3) == box .and. &
                   &abs(smpd_mask - smpd) <= 1.e-6
    end subroutine state_mask_is_compatible

    logical function should_regenerate_automask( params, which_iter, l_state_mask_exists, l_state_mask_compatible )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: which_iter
        logical,           intent(in) :: l_state_mask_exists, l_state_mask_compatible
        should_regenerate_automask = .false.
        if( trim(params%automsk).eq.'no' ) return
        if( .not. l_state_mask_exists ) then
            should_regenerate_automask = .true.
        else if( .not. l_state_mask_compatible ) then
            should_regenerate_automask = .true.
        else if( which_iter == params%startit ) then
            should_regenerate_automask = .true.
        else if( mod(which_iter, AMSK_FREQ) == 0 ) then
            should_regenerate_automask = .true.
        end if
    end function should_regenerate_automask

end module simple_vol_pproc_policy
