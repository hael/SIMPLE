!@descr: reusable construction of the abinitio3D docked multi-state split checkpoint
module simple_abinitio3D_split_checkpoint
use simple_commanders_api
use simple_abinitio_utils, only: abinitio_update_frac_max, calc_docked_multistate_max_sampling, &
    &calc_rec, cline_refine3D, randomize_states, set_cline_refine3D, update_frac
implicit none
#include "simple_local_flags.inc"

public :: build_abinitio3D_split_checkpoint
private

contains

    subroutine build_abinitio3D_split_checkpoint( params, spproj, xrefine3D, xrec3D, split_stage, &
        &nptcls_eff, nstates, force_full_sampling, post_split_update_frac )
        class(parameters),     intent(inout) :: params
        class(sp_project),     intent(inout) :: spproj
        class(commander_base), intent(inout) :: xrefine3D, xrec3D
        integer,               intent(in)    :: split_stage, nptcls_eff, nstates
        logical,               intent(in)    :: force_full_sampling
        real,                  intent(out)   :: post_split_update_frac
        integer :: cohort_target, nptcls_cap
        real    :: cohort_update_frac
        if( nptcls_eff < 1 ) THROW_HARD('abinitio3D split checkpoint requires active particles')
        if( nstates < 2 ) THROW_HARD('abinitio3D split checkpoint requires at least two states')
        if( force_full_sampling )then
            cohort_target      = nptcls_eff
            cohort_update_frac = 1.0
        else
            call calc_docked_multistate_max_sampling(params, nptcls_eff, nptcls_cap, cohort_update_frac)
            cohort_target = nptcls_cap
        endif
        call prepare_particle_cohort(params, spproj, xrefine3D, cohort_target, cohort_update_frac)
        params%nstates = nstates
        if( force_full_sampling )then
            post_split_update_frac = 1.0
        else
            post_split_update_frac = real(nstates * params%nsample) / real(nptcls_eff)
            post_split_update_frac = min(abinitio_update_frac_max(), post_split_update_frac)
        endif
        ! set_cline_refine3D reads the module-level update_frac, so it must hold
        ! the post-split value before the split-stage cline is emitted
        update_frac = post_split_update_frac
        call set_cline_refine3D(params, split_stage, l_cavgs=.false.)
        call randomize_states(params, spproj, params%projfile, xrec3D, split_stage, &
            &clean_sampling=.false., reconstruct_states=.false.)
        if( force_full_sampling )then
            call calc_rec(params, params%projfile, xrec3D, split_stage)
        else
            call select_reconstruction_sample(params, spproj, nptcls_eff, post_split_update_frac)
            call calc_rec(params, params%projfile, xrec3D, split_stage, current_sample_only=.true.)
        endif
        write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
            &'>>> ABINITIO3D SPLIT CHECKPOINT STAGE/NSTATES/POSTSPLIT_UPDATE_FRAC: ', &
            &split_stage, '/', params%nstates, '/', post_split_update_frac
    end subroutine build_abinitio3D_split_checkpoint

    subroutine prepare_particle_cohort( params, spproj, xrefine3D, nsample, ufrac )
        class(parameters),     intent(in)    :: params
        class(sp_project),     intent(inout) :: spproj
        class(commander_base), intent(inout) :: xrefine3D
        integer,               intent(in)    :: nsample
        real,                  intent(in)    :: ufrac
        integer, parameter :: COHORT_ASSIGNMENT_NITERS = 1
        type(cmdline) :: cline_assignment
        integer       :: iter_assignment, nactive, nupdated
        call spproj%read_segment('ptcl3D', params%projfile)
        call spproj%os_ptcl3D%clean_entry('updatecnt', 'sampled')
        call spproj%write_segment_inside(params%oritype, params%projfile)
        iter_assignment = next_refine3D_iteration()
        write(logfhandle,'(A,I0,A,I0,A,I0)') &
            &'>>> ABINITIO3D SPLIT COHORT TARGET/ACTIVE/ITER: ', &
            &nsample, '/', spproj%os_ptcl3D%count_state_gt_zero(), '/', iter_assignment
        cline_assignment = cline_refine3D
        call cline_assignment%set('prg',               'refine3D')
        call cline_assignment%set('mkdir',                   'no')
        call cline_assignment%set('refine',                'prob')
        call cline_assignment%set('balance',                'yes')
        call cline_assignment%set('nsample',              nsample)
        call cline_assignment%set('frac_best',                1.0)
        call cline_assignment%set('fillin',                  'no')
        call cline_assignment%set('update_frac',            ufrac)
        call cline_assignment%set('trail_rec',               'no')
        call cline_assignment%set('volrec',                  'no')
        call cline_assignment%set('sticky_class_sampling',   'no')
        call cline_assignment%set('maxits', COHORT_ASSIGNMENT_NITERS)
        call cline_assignment%set('startit',         iter_assignment)
        call cline_assignment%set('which_iter',      iter_assignment)
        call cline_assignment%set('extr_iter',       iter_assignment)
        call cline_assignment%delete('endit')
        call cline_assignment%delete('greedy_sampling')
        call xrefine3D%execute(cline_assignment)
        call cleanup_assignment_files(params%nparts)
        call read_assignment_coverage(params, spproj, nactive, nupdated)
        if( nupdated < nsample ) THROW_HARD('abinitio3D split cohort assignment updated too few particles')
        call cline_assignment%kill
    end subroutine prepare_particle_cohort

    subroutine select_reconstruction_sample( params, spproj, nptcls_eff, ufrac )
        class(parameters), intent(in)    :: params
        class(sp_project), intent(inout) :: spproj
        integer,           intent(in)    :: nptcls_eff
        real,              intent(in)    :: ufrac
        type(class_sample), allocatable :: split_clssmp(:)
        integer, allocatable :: split_pinds(:), sampled(:), states(:)
        integer :: noris, nrequested, nselected, sample_ind, state, state_pop
        call spproj%read_segment('ptcl3D', params%projfile)
        noris = spproj%os_ptcl3D%get_noris()
        call read_class_samples(split_clssmp, string(CLASS_SAMPLING_FILE))
        call spproj%os_ptcl3D%sample4update_class(split_clssmp, [1,noris], ufrac, &
            &nselected, split_pinds, .true., .false., sampled_only=.true.)
        nrequested = nint(ufrac * real(nptcls_eff))
        if( nselected /= nrequested )then
            THROW_HARD('abinitio3D split checkpoint failed to select the exact requested sample')
        endif
        call spproj%write_segment_inside(params%oritype, params%projfile)
        sampled   = spproj%os_ptcl3D%get_all_asint('sampled')
        states    = spproj%os_ptcl3D%get_all_asint('state')
        sample_ind = maxval(sampled)
        do state = 1, params%nstates
            state_pop = count(states == state .and. sampled == sample_ind)
            if( state_pop <= 5 )then
                THROW_HARD('abinitio3D split checkpoint sample has insufficient state population')
            endif
        enddo
        write(logfhandle,'(A,I0,A,I0,A,F8.4)') &
            &'>>> ABINITIO3D SPLIT RECONSTRUCTION SAMPLE SELECTED/COHORT/FRACTION: ', &
            &nselected, '/', count(sampled > 0), '/', real(nselected) / real(count(sampled > 0))
        if( allocated(sampled)     ) deallocate(sampled)
        if( allocated(states)      ) deallocate(states)
        if( allocated(split_pinds) ) deallocate(split_pinds)
        if( allocated(split_clssmp) ) call deallocate_class_samples(split_clssmp)
    end subroutine select_reconstruction_sample

    subroutine read_assignment_coverage( params, spproj, nactive, nupdated )
        class(parameters), intent(in)    :: params
        class(sp_project), intent(inout) :: spproj
        integer,           intent(out)   :: nactive, nupdated
        integer, allocatable :: states(:), updatecnts(:)
        call spproj%read_segment('ptcl3D', params%projfile)
        if( .not. spproj%os_ptcl3D%isthere('updatecnt') )then
            THROW_HARD('abinitio3D split checkpoint requires assignment update counts')
        endif
        states     = spproj%os_ptcl3D%get_all_asint('state')
        updatecnts = spproj%os_ptcl3D%get_all_asint('updatecnt')
        nactive    = count(states > 0)
        nupdated   = count(states > 0 .and. updatecnts > 0)
        if( allocated(states)     ) deallocate(states)
        if( allocated(updatecnts) ) deallocate(updatecnts)
    end subroutine read_assignment_coverage

    integer function next_refine3D_iteration() result(iter)
        iter = 1
        if( cline_refine3D%defined('endit') )then
            iter = cline_refine3D%get_iarg('endit') + 1
        else if( cline_refine3D%defined('which_iter') )then
            iter = cline_refine3D%get_iarg('which_iter') + 1
        endif
        iter = max(1, iter)
    end function next_refine3D_iteration

    subroutine cleanup_assignment_files( nparts )
        integer, intent(in) :: nparts
        call del_files(DIST_FBODY,       nparts, ext='.dat')
        call del_files(ASSIGNMENT_FBODY, nparts, ext='.dat')
        call del_file(DIST_FBODY//'.dat')
        call del_file(ASSIGNMENT_FBODY//'.dat')
    end subroutine cleanup_assignment_files

end module simple_abinitio3D_split_checkpoint
