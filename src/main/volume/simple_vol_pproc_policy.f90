!@descr: per-state automask policy decisions consumed by volume assembly
module simple_vol_pproc_policy
use simple_core_module_api
use simple_parameters, only: parameters
implicit none

private

public :: vol_pproc_plan
public :: plan_state_postprocess
public :: state_mask_is_compatible
public :: NU_ENVMASK_ACTION_NONE
public :: NU_ENVMASK_ACTION_REGENERATE

! Only two outcomes are acted on: regenerate this cycle, or keep whatever is on
! disk. There is no distinct "reuse" action because no consumer needs one.
integer, parameter :: NU_ENVMASK_ACTION_NONE       = 0
integer, parameter :: NU_ENVMASK_ACTION_REGENERATE = 1

type :: vol_pproc_plan
    type(string) :: nu_envmask_file
    integer      :: nu_envmask_action                  = NU_ENVMASK_ACTION_NONE
    logical      :: l_nu_envmask_enabled               = .false.
    logical      :: l_nu_envmask_exists                = .false.
    logical      :: l_nu_envmask_compatible            = .false.
    logical      :: l_nu_envmask_incompatible          = .false.
end type vol_pproc_plan

contains

    subroutine plan_state_postprocess( params, state, which_iter, plan )
        class(parameters),            intent(in)    :: params
        integer,                      intent(in)    :: state
        integer,                      intent(in)    :: which_iter
        type(vol_pproc_plan),         intent(inout) :: plan
        plan%nu_envmask_file = string(NU_ENVMASK_FBODY)//int2str_pad(state,2)//string(MRC_EXT)
        plan%nu_envmask_action                = NU_ENVMASK_ACTION_NONE
        plan%l_nu_envmask_enabled             = params%l_nonuniform .and. trim(params%automsk).ne.'no'
        plan%l_nu_envmask_exists              = .false.
        plan%l_nu_envmask_compatible          = .false.
        plan%l_nu_envmask_incompatible        = .false.
        if( .not. plan%l_nu_envmask_enabled ) return
        call state_mask_is_compatible(plan%nu_envmask_file, params%box_crop, params%smpd_crop, &
            &plan%l_nu_envmask_exists, plan%l_nu_envmask_compatible)
        plan%l_nu_envmask_incompatible = plan%l_nu_envmask_exists .and. (.not. plan%l_nu_envmask_compatible)
        if( should_regenerate_nu_envmask(params, which_iter, plan%l_nu_envmask_exists, &
            &plan%l_nu_envmask_compatible) )then
            plan%nu_envmask_action = NU_ENVMASK_ACTION_REGENERATE
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
        smpd_mask  = find_img_smpd(mskfile_state)
        compatible = all(ldim_mask == box) .and.  abs(smpd_mask - smpd) <= 1.e-6
    end subroutine state_mask_is_compatible

    logical function should_regenerate_nu_envmask( params, which_iter, l_state_mask_exists, l_state_mask_compatible )
        class(parameters), intent(in) :: params
        integer,           intent(in) :: which_iter
        logical,           intent(in) :: l_state_mask_exists, l_state_mask_compatible
        should_regenerate_nu_envmask = .false.
        if( trim(params%automsk).eq.'no' ) return
        if( .not. l_state_mask_exists ) then
            should_regenerate_nu_envmask = .true.
        else if( .not. l_state_mask_compatible ) then
            should_regenerate_nu_envmask = .true.
        else if( which_iter == params%startit ) then
            should_regenerate_nu_envmask = .true.
        else if( mod(which_iter, AMSK_FREQ) == 0 ) then
            should_regenerate_nu_envmask = .true.
        end if
    end function should_regenerate_nu_envmask

end module simple_vol_pproc_policy
