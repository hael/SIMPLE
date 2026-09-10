!@descr: validation and group loading of canonical sigma2 state.
!
!  Deliberately does NOT depend on simple_builder: builder depends on
!  simple_euclid_sigma2, so anything reachable from the sigma2 side must stay
!  builder-free. Callers pass the pieces (pftc, esig, oris) explicitly.
module simple_sigma2_files
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_euclid_sigma2, only: euclid_sigma2
use simple_parameters,    only: parameters
use simple_polarft_calc,  only: polarft_calc
use simple_sp_project,    only: sp_project
use simple_sigma2_state,  only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
use simple_sigma2_state_file, only: sigma2_state_validate_file, SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, &
    &SIGMA2_STATE_COMMITTED
implicit none

public :: canonical_sigma2_consumable, load_sigma2_groups
private
#include "simple_local_flags.inc"

contains


    !> Can a euclid consumer at the given native grid and grouping load the
    !! project's canonical sigma2 state: registered, file-valid, committed,
    !! same grid and shell range, same ordered particle layout, same grouping.
    !! The single definition shared by refine3D initialization, reconstruction
    !! loading, the final reconstruction and the sigma2 bootstrap service.
    logical function canonical_sigma2_consumable( project, os, box, smpd, l_sigma_glob, message ) result( l_ok )
        class(sp_project),     intent(inout) :: project
        class(oris),           intent(inout) :: os
        integer,               intent(in)    :: box
        real,                  intent(in)    :: smpd
        logical,               intent(in)    :: l_sigma_glob
        character(len=STDLEN), intent(out)   :: message
        type(string)   :: state_path
        integer(int64) :: layout_digest
        integer        :: expected_grouping, expected_ngroups, iptcl, noris, status
        logical        :: found
        l_ok    = .false.
        message = ''
        call project%get_sigma2_state_path(state_path, found)
        if( .not. found )then
            message = 'canonical sigma2 state is not registered in the project'
            return
        endif
        call sigma2_state_validate_file(state_path%to_char(), status, message, deep=.true.)
        if( status /= 0 )then
            call state_path%kill
            return
        endif
        noris         = os%get_noris()
        layout_digest = sigma2_state_project_layout_digest(project, os)
        if( layout_digest == 0_int64 )then
            message = 'canonical sigma2 layout digest is undefined for this project'
            call state_path%kill
            return
        endif
        if( l_sigma_glob )then
            expected_grouping = SIGMA2_GROUP_GLOBAL
            expected_ngroups  = 1
        else
            expected_grouping = SIGMA2_GROUP_STACK
            expected_ngroups  = 0
            do iptcl = 1, noris
                if( os%get_state(iptcl) <= 0 ) cycle
                expected_ngroups = max(expected_ngroups, os%get_int(iptcl, 'stkind'))
            enddo
        endif
        call sigma2_state_validate_identity(state_path%to_char(), box, smpd, 1, fdim(box)-1, noris, &
            &layout_digest, status, message, expected_state=SIGMA2_STATE_COMMITTED, &
            &expected_grouping=expected_grouping, expected_ngroups=expected_ngroups)
        l_ok = status == 0
        call state_path%kill
    end function canonical_sigma2_consumable

    !>  \brief  Loads grouped sigma2 from the validated project-registered state.
    !!
    !!          Only acts when objfun=euclid, mirroring the convention that
    !!          objfun=cc needs no sigmas.
    !!
    subroutine load_sigma2_groups( params, pftc, esig, project, os, loaded )
        class(parameters),    intent(inout) :: params
        class(polarft_calc),  intent(inout) :: pftc
        class(euclid_sigma2), intent(inout) :: esig
        class(sp_project),    intent(inout) :: project
        class(oris),          intent(inout) :: os
        logical,              intent(out)   :: loaded
        type(string) :: state_path
        integer      :: noris
        logical      :: state_path_found
        character(len=STDLEN) :: message
        loaded = .false.
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        noris = os%get_noris()
        if( noris < 1 ) return
        if( .not. canonical_sigma2_consumable(project, os, params%box, params%smpd, params%l_sigma_glob, message) ) &
            &THROW_HARD(trim(message))
        call project%get_sigma2_state_path(state_path, state_path_found)
        if( .not. state_path_found ) THROW_HARD('Canonical sigma2 state path disappeared after validation')
        call esig%new(params, pftc, state_path, params%box)
        call esig%read_groups(os)
        call esig%allocate_ptcls
        loaded = allocated(esig%sigma2_noise)
        call state_path%kill
    end subroutine load_sigma2_groups

end module simple_sigma2_files
