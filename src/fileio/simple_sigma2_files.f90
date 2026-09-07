!@descr: discovery and group loading of sigma2 files.
!
!  Workflows that consume sigma2 spectra estimated by an UPSTREAM run (rather
!  than estimating their own) pick the group STAR for the requested iteration
!  in the CURRENT execution directory and load per-particle spectra from it.
!  This module is the single owner of that logic.
!
!  There is deliberately NO implicit carry-over from other directories
!  (2026-09-06): the former sibling-directory discovery (glob over sibling
!  directories of the current run for the same project file)
!  copied every completed run of the same project into new runs and seeded
!  their later stage starts with foreign iteration-numbered STARs. refine3D
!  continue=yes copies its predecessor's partition files explicitly from the
!  recorded previous refinement directory; bootstrap_rec3D estimates its own
!  sigmas; a standalone euclid reconstruct3D must be run in, or pointed at,
!  the directory that holds the sigma files, and fails loudly otherwise.
!
!  Deliberately does NOT depend on simple_builder: builder depends on
!  simple_euclid_sigma2, so anything reachable from the sigma2 side must stay
!  builder-free. Callers pass the pieces (pftc, esig, oris) explicitly.
module simple_sigma2_files
use, intrinsic :: iso_fortran_env, only: int64
use simple_core_module_api
use simple_cmdline,       only: cmdline
use simple_euclid_sigma2, only: euclid_sigma2, sigma2_star_from_iter
use simple_parameters,    only: parameters
use simple_polarft_calc,  only: polarft_calc
use simple_sp_project,    only: sp_project
use simple_sigma2_state,  only: sigma2_state_project_layout_digest, sigma2_state_validate_identity
use simple_sigma2_state_file, only: SIGMA2_GROUP_GLOBAL, SIGMA2_GROUP_STACK, SIGMA2_STATE_COMMITTED
implicit none

public :: pick_sigma_group_file, load_sigma2_groups
private
#include "simple_local_flags.inc"

contains


    !>  \brief  Prefers the group star file for the current iteration, else the
    !!          highest-numbered sigma2 group star file present locally.
    subroutine pick_sigma_group_file( params, sigma_group_fname, found )
        class(parameters), intent(in)  :: params
        type(string),      intent(out) :: sigma_group_fname
        logical,           intent(out) :: found
        type(string), allocatable :: list(:)
        integer :: i, best_iter, iter_here
        sigma_group_fname = ''
        found             = .false.
        sigma_group_fname = sigma2_star_from_iter(params%which_iter)
        if( file_exists(sigma_group_fname) )then
            found = .true.
            return
        endif
        call sigma_group_fname%kill
        call simple_list_files(SIGMA2_GROUP_FBODY//'*'//STAR_EXT, list)
        if( .not. allocated(list) ) return
        if( size(list) < 1 )then
            deallocate(list)
            return
        endif
        best_iter = -1
        do i = 1, size(list)
            iter_here = trailing_iter_from_sigma_group_name(basename(list(i)))
            if( iter_here >= best_iter )then
                best_iter         = iter_here
                sigma_group_fname = basename(list(i))
                found             = .true.
            endif
        end do
        deallocate(list)
    end subroutine pick_sigma_group_file

    integer function trailing_iter_from_sigma_group_name( fname ) result( iter )
        type(string), intent(in) :: fname
        character(len=:), allocatable :: raw
        integer :: i, scale, d
        raw   = fname%to_char()
        iter  = 0
        scale = 1
        do i = len_trim(raw),1,-1
            if( raw(i:i) >= '0' .and. raw(i:i) <= '9' )then
                d     = iachar(raw(i:i)) - iachar('0')
                iter  = iter + d*scale
                scale = scale*10
            else if( scale > 1 )then
                exit
            endif
        end do
    end function trailing_iter_from_sigma_group_name

    !>  \brief  Loads grouped sigma2 into esig. Canonical mode resolves and
    !!          validates the project-registered committed state. Legacy mode
    !!          selects the group STAR in the current directory. Returns
    !!          whether usable spectra were obtained; a .false. result means
    !!          no STAR is present here -- callers decide what that means.
    !!
    !!          Only acts when objfun=euclid, mirroring the convention that
    !!          objfun=cc needs no sigmas.
    !!
    !!          LEGACY NOTE, because it is surprising: this temporarily widens
    !!          params%fromp/top to the whole particle range and forces global
    !!          sigma dispatch. Reused star files from upstream refine/abinitio
    !!          runs are usually single-group, so indexing groups by stkind
    !!          against a 1-group table would be wrong. fromp/top are restored
    !!          before returning; sigma_est is deliberately left set.
    subroutine load_sigma2_groups( params, pftc, esig, project, os, cline, loaded )
        class(parameters),    intent(inout) :: params
        class(polarft_calc),  intent(inout) :: pftc
        class(euclid_sigma2), intent(inout) :: esig
        class(sp_project),    intent(in)    :: project
        class(oris),          intent(inout) :: os
        class(cmdline),       intent(inout) :: cline
        logical,              intent(out)   :: loaded
        type(string)     :: sigma_part_fname, sigma_group_fname, state_path, cwd
        integer(int64)   :: layout_digest
        integer          :: fromp_saved, top_saved, noris, expected_grouping, expected_ngroups
        integer          :: iptcl, status
        logical          :: has_group, state_path_found
        character(len=STDLEN) :: message
        loaded = .false.
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        noris = os%get_noris()
        if( noris < 1 ) return
        if( params%l_sigma_canonical )then
            call project%get_sigma2_state_path(state_path, state_path_found)
            if( .not. state_path_found ) THROW_HARD('particle project has no canonical sigma2 state path')
            layout_digest = sigma2_state_project_layout_digest(project, os)
            if( params%l_sigma_glob )then
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
            call sigma2_state_validate_identity(state_path%to_char(), params%box, params%smpd, 1, &
                &fdim(params%box)-1, noris, layout_digest, status, message, &
                &expected_state=SIGMA2_STATE_COMMITTED, expected_grouping=expected_grouping, &
                &expected_ngroups=expected_ngroups)
            if( status /= 0 ) THROW_HARD(trim(message))
            call esig%new(params, pftc, state_path, params%box)
            call esig%read_groups(os)
            call esig%allocate_ptcls
            loaded = allocated(esig%sigma2_noise)
            call state_path%kill
            return
        endif
        call pick_sigma_group_file(params, sigma_group_fname, has_group)
        if( .not. has_group )then
            ! no implicit discovery elsewhere: the caller's directory must hold
            ! the sigma files (callers hard-error on loaded=.false.)
            call simple_getcwd(cwd)
            write(logfhandle,'(A)') '>>> SIGMA2: no grouped sigma STAR ('//SIGMA2_GROUP_FBODY//&
                &'<iter>'//STAR_EXT//') found in '//cwd%to_char()//&
                &'; euclid reconstruction must run in the directory that holds the sigma2 files'
            call cwd%kill
            return
        endif
        sigma_part_fname    = SIGMA2_FBODY//int2str_pad(params%part,params%numlen)//'.dat'
        fromp_saved         = params%fromp
        top_saved           = params%top
        params%sigma_est    = 'global'
        params%l_sigma_glob = .true.
        call cline%set('sigma_est','global')
        params%fromp = 1
        params%top   = noris
        call esig%new(params, pftc, sigma_part_fname, params%box)
        call esig%read_groups(os, fname=sigma_group_fname)
        call esig%allocate_ptcls()
        loaded       = allocated(esig%sigma2_noise)
        params%fromp = fromp_saved
        params%top   = top_saved
        call sigma_part_fname%kill
        call sigma_group_fname%kill
    end subroutine load_sigma2_groups

end module simple_sigma2_files
