!@descr: fixed-reference CC pose initialization shared by 3D refinement workflows
module simple_external_reference_pose_initialization
use simple_commanders_api
use simple_commanders_euclid, only: commander_calc_pspec, commander_calc_group_sigmas
use simple_estimate_ssnr,     only: lpstages_setlims
use simple_refine3D_fnames,   only: refine3D_startvol_fname, refine3D_state_vol_fname
implicit none
#include "simple_local_flags.inc"

public :: initialize_poses_against_external_references
private

contains

    subroutine initialize_poses_against_external_references( params, parent_cline, xrefine3D, xrec3D, &
        &nactive, reference_vols, checkpoint_vols, pose_init_iter )
        class(parameters),     intent(in)    :: params
        class(cmdline),        intent(inout) :: parent_cline
        class(commander_base), intent(inout) :: xrefine3D, xrec3D
        integer,               intent(in)    :: nactive
        type(string),          intent(in)    :: reference_vols(:)
        type(string),          intent(out)   :: checkpoint_vols(:)
        integer,               intent(in)    :: pose_init_iter
        integer, parameter :: NSAMPLE_POSE_INIT_CAP = 100000
        integer, parameter :: NSPACE_POSE_INIT       = 2500
        real,    parameter :: LP_POSE_INIT           = 15.0
        type(commander_calc_pspec)        :: xcalc_pspec
        type(commander_calc_group_sigmas) :: xcalc_group_sigmas
        type(cmdline) :: cline_sigma_bootstrap, cline_pose_init, cline_sigmas, cline_checkpoint
        type(lp_crop_inf) :: lpinfo_pose_init(1)
        type(string) :: startvol
        integer :: state, nsample_pose_init
        real    :: ufrac_pose_init
        if( nactive < 1 ) THROW_HARD('external-reference pose initialization requires active particles')
        if( size(reference_vols) < 1 ) THROW_HARD('external-reference pose initialization requires references')
        if( size(checkpoint_vols) /= size(reference_vols) )then
            THROW_HARD('external-reference pose-initialization checkpoint/reference state count mismatch')
        endif
        if( pose_init_iter < 1 ) THROW_HARD('external-reference pose-initialization iteration must be positive')
        nsample_pose_init = min(nactive, NSAMPLE_POSE_INIT_CAP)
        ufrac_pose_init   = real(nsample_pose_init) / real(nactive)
        call lpstages_setlims(params%box, 1, params%smpd, LP_POSE_INIT, LP_POSE_INIT, &
            &lpinfo_pose_init)
        cline_pose_init = parent_cline
        call cline_pose_init%set('prg',            'refine3D')
        call cline_pose_init%set('mkdir',                'no')
        call cline_pose_init%set('objfun',               'cc')
        call cline_pose_init%set('sigma_est',         'global')
        call cline_pose_init%set('refine',           'greedy')
        call cline_pose_init%set('greedy_sampling',     'yes')
        call cline_pose_init%set('balance',              'no')
        call cline_pose_init%set('trail_rec',            'no')
        call cline_pose_init%set('volrec',               'no')
        call cline_pose_init%set('cc_emit_sigma',        'yes')
        call cline_pose_init%set('update_frac', ufrac_pose_init)
        call cline_pose_init%set('nsample', nsample_pose_init)
        call cline_pose_init%set('nstates', size(reference_vols))
        call cline_pose_init%set('nspace', NSPACE_POSE_INIT)
        call cline_pose_init%set('maxits',                  1)
        call cline_pose_init%set('minits',                  1)
        call cline_pose_init%set('startit',  pose_init_iter)
        call cline_pose_init%set('which_iter', pose_init_iter)
        call cline_pose_init%set('extr_iter', pose_init_iter)
        call cline_pose_init%set('lp',          LP_POSE_INIT)
        call cline_pose_init%set('lpstop',      LP_POSE_INIT)
        call cline_pose_init%set('trs', lpinfo_pose_init(1)%trslim)
        do state = 1,size(reference_vols)
            if( .not. file_exists(reference_vols(state)) )then
                THROW_HARD('external-reference pose-initialization input volume does not exist')
            endif
            call cline_pose_init%set('vol'//int2str(state), reference_vols(state))
            startvol = refine3D_startvol_fname(state)
            if( trim(reference_vols(state)%to_char()) /= trim(startvol%to_char()) )then
                call simple_copy_file(reference_vols(state), startvol)
            endif
        enddo
        if( lpinfo_pose_init(1)%l_autoscale )then
            call cline_pose_init%set('box_crop',  lpinfo_pose_init(1)%box_crop)
            call cline_pose_init%set('smpd_crop', lpinfo_pose_init(1)%smpd_crop)
        else
            call cline_pose_init%delete('box_crop')
            call cline_pose_init%delete('smpd_crop')
        endif
        call prepare_pose_initialization_coverage(params%projfile)
        cline_sigma_bootstrap = parent_cline
        call cline_sigma_bootstrap%set('prg',                   'calc_pspec')
        call cline_sigma_bootstrap%set('mkdir',                         'no')
        call cline_sigma_bootstrap%set('objfun',                    'euclid')
        call cline_sigma_bootstrap%set('sigma_est',                 'global')
        call cline_sigma_bootstrap%set('cc_emit_sigma',                 'no')
        call cline_sigma_bootstrap%set('sigma_transition_ready',        'no')
        call cline_sigma_bootstrap%set('which_iter',          pose_init_iter)
        call cline_sigma_bootstrap%delete('part')
        call cline_sigma_bootstrap%delete('box_crop')
        call cline_sigma_bootstrap%delete('smpd_crop')
        call cline_sigma_bootstrap%delete('update_frac')
        call cline_sigma_bootstrap%delete('nsample')
        call cline_sigma_bootstrap%delete('endit')
        write(logfhandle,'(A,I0)') '>>> EXTERNAL-REFERENCE IMAGE-POWER SIGMA2 BOOTSTRAP ITERATION: ', pose_init_iter
        call xcalc_pspec%execute(cline_sigma_bootstrap)
        call cline_sigma_bootstrap%kill
        call cline_pose_init%delete('endit')
        write(logfhandle,'(A,I0,A,I0,A,F8.4,A,F6.1)') &
            &'>>> FIXED-REFERENCE CC POSE INITIALIZATION STATES/NSAMPLE/FRACTION/LP: ', &
            &size(reference_vols), '/', nsample_pose_init, '/', ufrac_pose_init, '/', LP_POSE_INIT
        call xrefine3D%execute(cline_pose_init)
        call validate_pose_initialized_states(params%projfile, size(reference_vols), nsample_pose_init)
        cline_sigmas = cline_pose_init
        call cline_sigmas%set('prg',        'calc_group_sigmas')
        call cline_sigmas%set('which_iter', pose_init_iter + 1)
        call xcalc_group_sigmas%execute(cline_sigmas)
        call cline_sigmas%kill
        cline_checkpoint = cline_pose_init
        call cline_checkpoint%set('prg',         'reconstruct3D')
        call cline_checkpoint%set('objfun',      'cc')
        call cline_checkpoint%set('postprocess', 'no')
        call cline_checkpoint%set('nu_refine',   'no')
        call cline_checkpoint%delete('trail_rec')
        call cline_checkpoint%delete('refine')
        call cline_checkpoint%delete('objfun_den')
        call cline_checkpoint%delete('objfun_den_w')
        call cline_checkpoint%delete('sigma_est')
        call cline_checkpoint%delete('cc_emit_sigma')
        call cline_checkpoint%delete('ufrac_trec')
        call cline_checkpoint%delete('box_crop')
        call cline_checkpoint%delete('smpd_crop')
        call cline_checkpoint%delete('endit')
        call xrec3D%execute(cline_checkpoint)
        do state = 1,size(reference_vols)
            checkpoint_vols(state) = refine3D_state_vol_fname(state)
        enddo
        call cleanup_assignment_files(params%nparts)
        call cline_checkpoint%kill
        call cline_pose_init%kill
        call startvol%kill
        call parent_cline%set('sigma_transition_ready', 'yes')
        write(logfhandle,'(A)') &
            &'>>> FIXED-REFERENCE CC POSE INITIALIZATION COMPLETE; DATA-DERIVED CHECKPOINT MAPS ARE AUTHORITATIVE'
    end subroutine initialize_poses_against_external_references

    subroutine prepare_pose_initialization_coverage( projfile )
        type(string), intent(in) :: projfile
        type(sp_project) :: pose_proj
        call pose_proj%read_segment('ptcl3D', projfile)
        call pose_proj%os_ptcl3D%clean_entry('sampled', 'updatecnt')
        call pose_proj%write_segment_inside('ptcl3D', projfile)
        call pose_proj%kill
    end subroutine prepare_pose_initialization_coverage

    subroutine validate_pose_initialized_states( projfile, nstates, nsample_expected )
        type(string), intent(in) :: projfile
        integer,      intent(in) :: nstates, nsample_expected
        type(sp_project) :: pose_proj
        integer, allocatable :: states(:), updatecnts(:)
        integer :: state, nupdated, state_pop, nexpected
        call pose_proj%read_segment('ptcl3D', projfile)
        if( .not. pose_proj%os_ptcl3D%isthere('updatecnt') )then
            call pose_proj%kill
            THROW_HARD('external-reference pose initialization did not record update coverage')
        endif
        states     = pose_proj%os_ptcl3D%get_all_asint('state')
        updatecnts = pose_proj%os_ptcl3D%get_all_asint('updatecnt')
        nupdated   = count(states > 0 .and. updatecnts > 0)
        nexpected  = floor(0.99*real(nsample_expected)) ! to account for rouding
        if( nupdated < nexpected )then
            THROW_HARD('external-reference pose initialization updated fewer particles than the capped cohort')
        endif
        do state = 1,nstates
            state_pop = count(states == state .and. updatecnts > 0)
            write(logfhandle,'(A,I0,A,I0)') &
                &'>>> FIXED-REFERENCE CC POSE-INITIALIZED STATE/COHORT POPULATION: ', state, '/', state_pop
            if( state_pop < 1 )then
                THROW_HARD('external-reference pose initialization produced an empty state')
            endif
        enddo
        if( allocated(states)     ) deallocate(states)
        if( allocated(updatecnts) ) deallocate(updatecnts)
        call pose_proj%kill
    end subroutine validate_pose_initialized_states

    subroutine cleanup_assignment_files( nparts )
        integer, intent(in) :: nparts
        call del_files(DIST_FBODY,       nparts, ext='.dat')
        call del_files(ASSIGNMENT_FBODY, nparts, ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
    end subroutine cleanup_assignment_files

end module simple_external_reference_pose_initialization
