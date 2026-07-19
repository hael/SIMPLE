module simple_commanders_prob
use simple_commanders_api
use simple_pftc_srch_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_prob_tab
  contains
    procedure :: execute      => exec_prob_tab
end type commander_prob_tab

type, extends(commander_base) :: commander_prob_tab_neigh
    contains
        procedure :: execute      => exec_prob_tab_neigh
end type commander_prob_tab_neigh

type, extends(commander_base) :: commander_prob_align
  contains
    procedure :: execute      => exec_prob_align
end type commander_prob_align

type, extends(commander_base) :: commander_prob_align_neigh
    contains
        procedure :: execute      => exec_prob_align_neigh
end type commander_prob_align_neigh

type, extends(commander_base) :: commander_prob_tab2D
  contains
    procedure :: execute      => exec_prob_tab2D
end type commander_prob_tab2D

type, extends(commander_base) :: commander_prob_align2D
  contains
    procedure :: execute      => exec_prob_align2D
end type commander_prob_align2D

contains

    subroutine exec_prob_tab( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab,            only: eul_prob_tab
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab)       :: eulprob_obj_part
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz
        integer, allocatable :: batches(:,:)
        logical :: l_state_only
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,need3Dobjs=.false.)
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        if( nptcls < 1 ) THROW_HARD('exec_prob_tab selected no particles')
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call read_reprojection_model(params, build, batchsz_max)
        call prep_sigmas_objfun(params, build, .false.)
        call alloc_ptcl_imgs( params, build, tmp_imgs, tmp_imgs_pad, batchsz_max )
        call build%pftc%memoize_refs(eulspace=build%eulspace)
        ! Fill the partition table in matcher-sized batches to cap particle PFT memo memory.
        l_state_only = str_has_substr(params%refine, 'prob_state')
        call eulprob_obj_part%new_worker(params,build,pinds)
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%begin_write(fname)
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad)
            if( l_state_only )then
                call eulprob_obj_part%fill_tab_state_only_range(batch_start, batch_end)
            else
                call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
            endif
        end do
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call build%pftc%kill
        call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_matcher_2Dprep
        use simple_matcher_refvol_utils,    only: read_reprojection_model
        use simple_matcher_ptcl_batch,      only: prep_sigmas_objfun, alloc_ptcl_imgs, build_batch_particles3D, clean_batch_particles3D
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(image), allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab_neigh) :: eulprob_obj_part_neigh
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,need3Dobjs=.true.)
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        if( nptcls < 1 ) THROW_HARD('exec_prob_tab_neigh selected no particles')
        ! All neighborhood modes can fill the table in matcher-sized batches; this
        ! caps particle-image/PFT memo memory without changing assignment ownership.
        fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        call run_prob_tab_neigh_batch(fname)
        call fname%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine prepare_prob_neigh_workspace(batchsz_here)
            integer, intent(in) :: batchsz_here
            call read_reprojection_model(params, build, batchsz_here)
            call prep_sigmas_objfun(params, build, .false.)
            call alloc_ptcl_imgs(params, build, tmp_imgs, tmp_imgs_pad, batchsz_here)
            call build%pftc%memoize_refs(eulspace=build%eulspace)
        end subroutine prepare_prob_neigh_workspace

        subroutine cleanup_prob_neigh_workspace
            call build%pftc%kill
            call clean_batch_particles3D(build, tmp_imgs, tmp_imgs_pad)
        end subroutine cleanup_prob_neigh_workspace

        subroutine run_prob_tab_neigh_batch(outfname)
            class(string), intent(in)  :: outfname
            integer, allocatable :: batches(:,:)
            integer :: ibatch, batch_start, batch_end, batchsz, batchsz_max, nbatches
            batchsz_max = min(nptcls, max(1, params%nthr * BATCHTHRSZ))
            nbatches    = ceiling(real(nptcls) / real(batchsz_max))
            batches     = split_nobjs_even(nptcls, nbatches)
            batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
            call prepare_prob_neigh_workspace(batchsz_max)
            call eulprob_obj_part_neigh%new_neigh(params,build,pinds)
            call eulprob_obj_part_neigh%begin_write(outfname)
            do ibatch = 1, nbatches
                batch_start = batches(ibatch,1)
                batch_end   = batches(ibatch,2)
                batchsz     = batch_end - batch_start + 1
                call build_batch_particles3D(params, build, batchsz, pinds(batch_start:batch_end), tmp_imgs, tmp_imgs_pad)
                call eulprob_obj_part_neigh%fill_tab_range(batch_start, batch_end)
            enddo
            call eulprob_obj_part_neigh%write_tab(outfname)
            call eulprob_obj_part_neigh%kill
            if( allocated(batches) ) deallocate(batches)
            call cleanup_prob_neigh_workspace
        end subroutine run_prob_tab_neigh_batch

    end subroutine exec_prob_tab_neigh


    subroutine exec_prob_align( self, cline )
        use simple_eul_prob_tab,            only: eul_prob_tab
        use simple_prob_posterior3D,          only: POSTERIOR3D_FNAME
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(commander_prob_tab) :: xprob_tab
        type(eul_prob_tab)       :: eulprob_obj_glob
        type(cmdline)            :: cline_prob_tab
        type(qsys_env)           :: qenv
        type(chash)              :: job_descr
        integer :: nptcls, ipart
        logical :: l_state_only
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, need3Dobjs=.false.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .false.)
        ! generating all corrs on all parts
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab' ) ! required for distributed call
        ! execution
        if( .not.cline_prob_tab%defined('nparts') )then
            call xprob_tab%execute(cline_prob_tab)
        else
            ! setup the environment for distributed execution
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            ! schedule
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! Build the global table only after worker tables are complete.  Keeping it
        ! live while workers build dense partition tables roughly doubles peak RSS.
        l_state_only = str_has_substr(params%refine, 'prob_state')
        if( l_state_only )then
            call eulprob_obj_glob%new_state(params,build,pinds)
        else
            call eulprob_obj_glob%new(params,build,pinds)
        endif
        ! reading corrs from all parts
        if( l_state_only )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_state_tab(fname)
            enddo
            call eulprob_obj_glob%state_assign
        else
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
                call eulprob_obj_glob%read_tab_to_glob(fname)
            enddo
            call eulprob_obj_glob%ref_assign
        endif
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob%write_assignment(fname)
        if( trim(params%write_prior) == 'yes' )then
            ! The dense worker path intentionally avoids allocating 3D volume
            ! objects.  The artifact only needs the lightweight projection
            ! orientation grid for source Euler coordinates.
            call build%eulspace%new(params%nspace, is_ptcl=.false.)
            call build%pgrpsyms%build_refspiral(build%eulspace)
            call eulprob_obj_glob%write_posterior3D(string(POSTERIOR3D_FNAME))
        endif
        call eulprob_obj_glob%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

    subroutine exec_prob_align_neigh( self, cline )
        use simple_eul_prob_tab_neigh,      only: eul_prob_tab_neigh
        use simple_prob_posterior3D,        only: POSTERIOR3D_FNAME
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4fillin, sample_ptcls4update3D
        use simple_builder,                 only: builder
        class(commander_prob_align_neigh), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        integer,           allocatable :: pinds(:)
        type(string)                   :: fname
        type(builder)                  :: build
        type(parameters)               :: params
        type(commander_prob_tab_neigh) :: xprob_tab_neigh
        type(eul_prob_tab_neigh)       :: eulprob_obj_glob_neigh
        type(cmdline)                  :: cline_prob_tab
        type(qsys_env)                 :: qenv
        type(chash)                    :: job_descr
        integer :: nptcls
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, need3Dobjs=.true.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update3D(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        call build%spproj%write_segment_inside(params%oritype)
        call cleanup_prob_align_outputs(params, .true.)
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab_neigh')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab_neigh%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        ! Construct global storage only after worker scoring has released its table.
        call eulprob_obj_glob_neigh%new_neigh_global(params,build,pinds)
        call eulprob_obj_glob_neigh%read_tabs_to_glob(string(DIST_FBODY)//'_neigh_', params%nparts, params%numlen)
        call eulprob_obj_glob_neigh%ref_assign
        if( trim(params%prob_neigh_mode) == 'posterior' )then
            ! The merged current candidate table is the evidence for the next
            ! posterior generation.  Publish only after global assignment has
            ! consumed the table for the current iteration.
            call eulprob_obj_glob_neigh%write_posterior3D_candidates(POSTERIOR3D_FNAME)
        endif
        ! write the iptcl->(iref,istate) assignment
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        call eulprob_obj_glob_neigh%write_assignment(fname)
        call eulprob_obj_glob_neigh%kill
        ! cleanup
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align_neigh'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN_NEIGH NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align_neigh

    subroutine exec_prob_tab2D( self, cline )
        use simple_matcher_smpl_and_lplims, only: set_bp_range2D
        use simple_strategy2D_matcher,  only: set_b_p_ptrs2D, &
                                              ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad
        use simple_matcher_pftc_prep,      only: prep_pftc4align2D
        use simple_matcher_ptcl_batch,  only: alloc_ptcl_imgs, build_batch_particles2D, clean_batch_particles2D
        use simple_imgarr_utils,        only: alloc_imgarr
        use simple_classaverager,       only: cavger_new, cavger_read_all, cavger_kill
        use simple_eul_prob_tab2D,      only: eul_prob_tab2D
        class(commander_prob_tab2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,     allocatable :: pinds(:)
        type(string)             :: fname
        type(builder)            :: build
        type(parameters)         :: params
        type(eul_prob_tab2D)     :: eulprob_obj_part
        real    :: frac_srch_space
        integer :: nptcls, batchsz_max, nbatches, ibatch, batch_start, batch_end, batchsz
        integer, allocatable :: batches(:,:)
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline, params, need3Dobjs=.false.)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        frac_srch_space  = build%spproj_field%get_avg('frac')
        call set_bp_range2D(params, build, cline, params%which_iter, frac_srch_space)
        ! reproduce particle sampling from exec_prob_align2D
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab2D requires prior particle sampling (in exec_prob_align2D)')
        endif
        batchsz_max = min(nptcls, params%nthr * BATCHTHRSZ)
        nbatches    = ceiling(real(nptcls) / real(batchsz_max))
        batches     = split_nobjs_even(nptcls, nbatches)
        batchsz_max = maxval(batches(:,2) - batches(:,1) + 1)
        call set_b_p_ptrs2D(params, build)
        call alloc_ptcl_imgs(params, build, ptcl_match_imgs, ptcl_match_imgs_pad, batchsz_max)
        call alloc_imgarr(batchsz_max, [params%box, params%box, 1], params%smpd, ptcl_imgs)
        ! mirror cluster2D_exec reference setup
        call cavger_new(params, build)
        if( .not. cline%defined('refs') ) THROW_HARD('exec_prob_tab2D requires refs on the command line')
        call cavger_read_all
        call prep_pftc4align2D(params, build, ptcl_match_imgs_pad, batchsz_max, params%which_iter, .false.)
        ! Fill the partition table in matcher-sized batches to cap polar FT memo memory.
        call eulprob_obj_part%new_worker(params,build,pinds)
        do ibatch = 1, nbatches
            batch_start = batches(ibatch,1)
            batch_end   = batches(ibatch,2)
            batchsz     = batch_end - batch_start + 1
            call build_batch_particles2D(params, build, batchsz, pinds(batch_start:batch_end),&
                &ptcl_imgs(1:batchsz), ptcl_match_imgs, ptcl_match_imgs_pad)
            call eulprob_obj_part%fill_tab_range(batch_start, batch_end)
        end do
        ! write the 2D probability table
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part%write_tab(fname)
        call eulprob_obj_part%kill
        if( allocated(batches) ) deallocate(batches)
        call clean_batch_particles2D(build, ptcl_imgs, ptcl_match_imgs, ptcl_match_imgs_pad)
        call cavger_kill
        call build%pftc%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_tab2D'))
        call simple_end('**** SIMPLE_PROB_TAB2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab2D

    subroutine exec_prob_align2D( self, cline )
        use simple_eul_prob_tab2D,          only: eul_prob_tab2D, PRIOR2D_STAGE5_FNAME
        use simple_strategy2D_matcher,      only: set_b_p_ptrs2D
        use simple_matcher_smpl_and_lplims, only: sample_ptcls4update2D
        use simple_builder,                 only: builder
        class(commander_prob_align2D), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        integer,       allocatable :: pinds(:)
        type(string)               :: fname
        type(builder)              :: build
        type(parameters)           :: params
        type(commander_prob_tab2D) :: xprob_tab2D
        type(eul_prob_tab2D)       :: eulprob_obj_glob
        type(cmdline)              :: cline_prob_tab
        type(qsys_env)             :: qenv
        type(chash)                :: job_descr
        integer :: nptcls, ipart
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, need3Dobjs=.false.)
        call set_b_p_ptrs2D(params, build)
        if( build%spproj_field%get_nevenodd() == 0 )then
            call build%spproj_field%partition_eo
            call build%spproj%write_segment_inside(params%oritype, params%projfile)
        endif
        if( params%startit == 1 .and. params%which_iter == params%startit )then
            call build%spproj_field%clean_entry('updatecnt', 'sampled')
        endif
        ! Mirror the 3D workflow: sampled-update is active from the first stage onward.
        ! In probabilistic mode the sampled subset is reused within the current iteration
        ! by prob_tab2D/cluster2D_exec, but it is redrawn on later iterations.
        call sample_ptcls4update2D(params, build, [params%fromp,params%top], params%l_update_frac, nptcls, pinds)
        write(logfhandle,'(A,I0,A,I0,A,I0)') '>>> PROB_ALIGN2D: sampled ', nptcls, ' particles over ', params%nparts, ' part(s)'
        call flush(logfhandle)
        ! write sampling to project
        call build%spproj%write_segment_inside(params%oritype)
        ! build the global prob table (nclasses x nptcls)
        call eulprob_obj_glob%new(params, build, pinds)
        ! generate partition-wise dist tables
        cline_prob_tab = cline
        call cline_prob_tab%set('prg', 'prob_tab2D')
        if( .not. cline_prob_tab%defined('nparts') )then
            call xprob_tab2D%execute(cline_prob_tab)
        else
            call qenv%new(params, params%nparts, nptcls=params%nptcls)
            call cline_prob_tab%gen_job_descr(job_descr)
            call qenv%gen_scripts_and_schedule_jobs(job_descr, array=L_USE_SLURM_ARR, extra_params=params)
        endif
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: prob_tab2D workers completed; merging partition tables'
        call flush(logfhandle)
        ! merge all partition tables into global
        do ipart = 1, params%nparts
            fname = string(DIST_FBODY)//int2str_pad(ipart,params%numlen)//'.dat'
            call eulprob_obj_glob%read_tab_to_glob(fname)
        end do
        ! global probabilistic class assignment
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: running global probabilistic assignment'
        call flush(logfhandle)
        call eulprob_obj_glob%ref_assign
        ! write assignment to file
        fname = string(ASSIGNMENT_FBODY)//'.dat'
        write(logfhandle,'(A,A)') '>>> PROB_ALIGN2D: writing assignment ', fname%to_char()
        call flush(logfhandle)
        call eulprob_obj_glob%write_assignment(fname)
        write(logfhandle,'(A)') '>>> PROB_ALIGN2D: assignment written'
        call flush(logfhandle)
        ! write per-particle prior ranking only when the controller has flagged this as the
        ! prior-production stage (stage PROB_PRIOR_STAGE-1, i.e. stage 5 by default)
        if( trim(params%write_prior) == 'yes' )then
            fname = string(PRIOR2D_STAGE5_FNAME)
            call eulprob_obj_glob%write_prior_topk(fname)
            write(logfhandle,'(A,A)') '>>> PROB_ALIGN2D: prior ranking written ', fname%to_char()
            call flush(logfhandle)
        endif
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_prob :: exec_prob_align2D'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN2D NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align2D

    subroutine cleanup_prob_align_outputs( params, neigh )
        class(parameters), intent(in) :: params
        logical,           intent(in) :: neigh
        type(string) :: fname
        integer :: ipart
        if( neigh )then
            do ipart = 1, params%nparts
                fname = string(DIST_FBODY)//'_neigh_'//int2str_pad(ipart,params%numlen)//'.dat'
                call del_file(fname)
            end do
        else
            call del_files(DIST_FBODY, params%nparts, ext='.dat', numlen=params%numlen)
        endif
        call del_file(string(ASSIGNMENT_FBODY)//'.dat')
        call fname%kill
    end subroutine cleanup_prob_align_outputs

end module simple_commanders_prob
