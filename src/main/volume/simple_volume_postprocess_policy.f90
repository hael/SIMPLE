module simple_volume_postprocess_policy
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

public :: volume_postprocess_plan
public :: plan_state_postprocess
public :: state_mask_is_compatible
private

type :: volume_postprocess_plan
    type(string) :: mskfile_state
    logical      :: l_automask_enabled            = .false.
    logical      :: l_state_mask_exists           = .false.
    logical      :: l_state_mask_compatible       = .false.
    logical      :: l_state_mask_incompatible     = .false.
    logical      :: regenerate_automask           = .false.
    logical      :: use_state_mask_for_nonuniform = .false.
end type volume_postprocess_plan

contains

    subroutine plan_state_postprocess( params, state, which_iter, l_nonuniform_mode, plan )
        class(parameters),            intent(in)    :: params
        integer,                      intent(in)    :: state
        integer,                      intent(in)    :: which_iter
        logical,                      intent(in)    :: l_nonuniform_mode
        type(volume_postprocess_plan), intent(inout) :: plan
        plan%mskfile_state             = string(AUTOMASK_FBODY)//int2str_pad(state,2)//string(MRC_EXT)
        plan%l_automask_enabled        = trim(params%automsk).ne.'no'
        plan%l_state_mask_exists       = .false.
        plan%l_state_mask_compatible   = .false.
        plan%l_state_mask_incompatible = .false.
        plan%regenerate_automask       = .false.
        plan%use_state_mask_for_nonuniform = .false.
        if( .not. plan%l_automask_enabled ) return
        call state_mask_is_compatible(plan%mskfile_state, params%box_crop, params%smpd_crop, &
            &plan%l_state_mask_exists, plan%l_state_mask_compatible)
        plan%l_state_mask_incompatible = plan%l_state_mask_exists .and. (.not. plan%l_state_mask_compatible)
        plan%regenerate_automask       = should_regenerate_automask(params, which_iter, &
                                        &plan%l_state_mask_exists, plan%l_state_mask_compatible)
        plan%use_state_mask_for_nonuniform = l_nonuniform_mode .and. &
            &(plan%regenerate_automask .or. plan%l_state_mask_compatible)
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
        call find_ldim_nptcls(mskfile_state, ldim_mask, nptcls_mask, smpd=smpd_mask)
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

end module simple_volume_postprocess_policy
