!@descr: discovery, carry-over and group loading of sigma2 files.
!
!  Workflows that consume sigma2 spectra estimated by an UPSTREAM run (rather
!  than estimating their own) all need the same three things: pull any sigma2
!  files left by a previous refine/abinitio run into the current execution
!  directory, pick the most recent group star file, and load per-particle
!  spectra from it. That logic originated inside simple_flex_analysis_strategy
!  and was subsequently duplicated into the experimental PCG reconstruction
!  commander; this module is the single owner so it cannot drift again.
!
!  Deliberately does NOT depend on simple_builder: builder depends on
!  simple_euclid_sigma2, so anything reachable from the sigma2 side must stay
!  builder-free. Callers pass the pieces (pftc, esig, oris) explicitly.
module simple_sigma2_files
use simple_core_module_api
use simple_cmdline,       only: cmdline
use simple_euclid_sigma2, only: euclid_sigma2, sigma2_star_from_iter
use simple_parameters,    only: parameters
use simple_polarft_calc,  only: polarft_calc
use simple_sp_project,    only: sp_project
implicit none

public :: carry_over_prior_sigma_files, pick_sigma_group_file, load_sigma2_groups
private
#include "simple_local_flags.inc"

contains

    !>  \brief  Copies sigma2 files left by prior runs into the current working
    !!          directory, searching sibling run directories first (e.g.
    !!          ../5_abinitio3D/<projfile>), then the project's recorded cwd,
    !!          the project's own directory, and the parent. Never overwrites a
    !!          file that is already local.
    subroutine carry_over_prior_sigma_files( params )
        class(parameters), intent(in) :: params
        type(sp_project)          :: inproj
        type(string), allocatable :: siblings(:)
        type(string) :: proj_dir, parent_dir, proj_cwd, proj_base, sibling_dir
        logical      :: have_proj_cwd
        integer      :: i
        proj_dir      = get_fpath(params%projfile)
        parent_dir    = string(PATH_PARENT)
        proj_base     = basename(params%projfile)
        proj_cwd      = ''
        have_proj_cwd = .false.
        ! Prefer sibling-run discovery from ../, e.g. ../5_abinitio3D/<projfile>.
        if( parent_dir /= '' .and. proj_base /= '' )then
            call simple_list_files(parent_dir%to_char()//'*/'//proj_base%to_char(), siblings)
            if( allocated(siblings) )then
                do i = 1, size(siblings)
                    sibling_dir = get_fpath(siblings(i))
                    if( sibling_dir /= '' ) call copy_sigma_files_from_dir(sibling_dir)
                    call sibling_dir%kill
                end do
                deallocate(siblings)
            endif
        endif
        if( file_exists(params%projfile) )then
            call inproj%read_non_data_segments(params%projfile)
            if( inproj%projinfo%isthere('cwd') )then
                call inproj%projinfo%getter(1,'cwd',proj_cwd)
                have_proj_cwd = proj_cwd /= ''
            endif
            call inproj%kill
        endif
        if( have_proj_cwd    ) call copy_sigma_files_from_dir(proj_cwd)
        if( proj_dir   /= '' ) call copy_sigma_files_from_dir(proj_dir)
        if( parent_dir /= '' ) call copy_sigma_files_from_dir(parent_dir)
        call proj_base%kill
        call proj_cwd%kill
        call proj_dir%kill
        call parent_dir%kill
    end subroutine carry_over_prior_sigma_files

    subroutine copy_sigma_files_from_dir( src_dir )
        type(string), intent(in) :: src_dir
        type(string), allocatable     :: list(:)
        type(string)                  :: target_name
        character(len=:), allocatable :: src_prefix
        integer :: i
        src_prefix = trim(src_dir%to_char())
        if( len_trim(src_prefix) < 1 ) return
        if( src_prefix(len_trim(src_prefix):len_trim(src_prefix)) /= '/' ) src_prefix = src_prefix//'/'
        call simple_list_files(src_prefix//SIGMA2_FBODY//'*', list)
        if( allocated(list) )then
            do i = 1, size(list)
                target_name = string(PATH_HERE)//basename(list(i))
                if( .not. file_exists(target_name) ) call simple_copy_file(list(i),target_name)
                call target_name%kill
            end do
            deallocate(list)
        endif
        call simple_list_files(src_prefix//SIGMA2_GROUP_FBODY//'*'//STAR_EXT, list)
        if( allocated(list) )then
            do i = 1, size(list)
                target_name = string(PATH_HERE)//basename(list(i))
                if( .not. file_exists(target_name) ) call simple_copy_file(list(i),target_name)
                call target_name%kill
            end do
            deallocate(list)
        endif
    end subroutine copy_sigma_files_from_dir

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

    !>  \brief  Carries over prior sigma2 files, picks the group star file and
    !!          loads per-particle-per-shell spectra into esig. Returns whether
    !!          usable spectra were obtained; a .false. result is a normal
    !!          outcome (no upstream sigma2 available), not an error -- callers
    !!          decide what that means for them.
    !!
    !!          Only acts when objfun=euclid, mirroring the convention that
    !!          objfun=cc needs no sigmas.
    !!
    !!          NOTE, because it is surprising: this temporarily widens
    !!          params%fromp/top to the whole particle range and forces global
    !!          sigma dispatch. Reused star files from upstream refine/abinitio
    !!          runs are usually single-group, so indexing groups by stkind
    !!          against a 1-group table would be wrong. fromp/top are restored
    !!          before returning; sigma_est is deliberately left set.
    subroutine load_sigma2_groups( params, pftc, esig, os, cline, loaded )
        class(parameters),    intent(inout) :: params
        class(polarft_calc),  intent(inout) :: pftc
        class(euclid_sigma2), intent(inout) :: esig
        class(oris),          intent(inout) :: os
        class(cmdline),       intent(inout) :: cline
        logical,              intent(out)   :: loaded
        type(string) :: sigma_part_fname, sigma_group_fname
        integer      :: fromp_saved, top_saved, noris
        logical      :: has_group
        loaded = .false.
        call carry_over_prior_sigma_files(params)
        call pick_sigma_group_file(params, sigma_group_fname, has_group)
        if( .not. has_group ) return
        if( params%cc_objfun /= OBJFUN_EUCLID ) return
        noris = os%get_noris()
        if( noris < 1 ) return
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
