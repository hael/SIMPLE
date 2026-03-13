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

contains

    subroutine exec_prob_tab( self, cline )
        use simple_strategy2D3D_common
        use simple_eul_prob_tab, only: eul_prob_tab
        use simple_imgarr_utils, only: dealloc_imgarr
        class(commander_prob_tab), intent(inout) :: self
        class(cmdline),            intent(inout) :: cline
        integer,          allocatable :: pinds(:)
        type(image),      allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)                  :: fname
        type(builder)                 :: build
        type(parameters)              :: params
        type(eul_prob_tab)            :: eulprob_obj_part
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( params, build, cline )
        ! The policy here ought to be that nothing is done with regards to sampling other than reproducing
        ! what was generated in the driver (prob_align, below). Sampling is delegated to prob_align (below)
        ! and merely reproduced here
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab requires prior particle sampling (in exec_prob_align)')
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=(params%l_polar .and. (.not.cline%defined('vol1'))) )           
        ! Build polar particle images
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        ! Filling prob table in eul_prob_tab
        call eulprob_obj_part%new(params, build, pinds)
        fname = string(DIST_FBODY)//int2str_pad(params%part,params%numlen)//'.dat'
        if( str_has_substr(params%refine, 'prob_state') )then
            call eulprob_obj_part%fill_tab_state_only
            call eulprob_obj_part%write_state_tab(fname)
        else
            call eulprob_obj_part%fill_tab
            call eulprob_obj_part%write_tab(fname)
        endif
        call eulprob_obj_part%kill
        call killimgbatch(build)
        call build%pftc%kill
        call build%kill_general_tbox
        call dealloc_imgarr(tmp_imgs)
        call dealloc_imgarr(tmp_imgs_pad)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab'))
        call simple_end('**** SIMPLE_PROB_TAB NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_tab

    subroutine exec_prob_tab_neigh( self, cline )
        use simple_strategy2D3D_common
        use simple_eul_prob_tab_neigh, only: eul_prob_tab_neigh
        use simple_eul_prob_tab,       only: angle_sampling, calc_athres, eulprob_dist_switch
        use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
        use simple_eulspace_neigh_map, only: eulspace_neigh_map
        use simple_imgarr_utils,       only: dealloc_imgarr
        class(commander_prob_tab_neigh), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer,            allocatable :: pinds(:)
        integer,            allocatable :: ptcl_state(:)
        logical,            allocatable :: neigh_mask(:,:)
        type(image),        allocatable :: tmp_imgs(:), tmp_imgs_pad(:)
        type(string)                    :: fname
        type(builder)                   :: build
        type(parameters)                :: params
        type(eul_prob_tab_neigh)        :: eulprob_obj_part_neigh
        type(eulspace_neigh_map)        :: neigh_map
        integer :: nptcls
        call cline%set('mkdir', 'no')
        call build%init_params_and_build_general_tbox(cline,params,do3d=.true.)
        call set_bp_range( params, build, cline )
        ! Sampling policy mirrors exec_prob_tab: only reproduce already sampled particles.
        if( build%spproj_field%has_been_sampled() )then
            call build%spproj_field%sample4update_reprod([params%fromp,params%top], nptcls, pinds)
        else
            THROW_HARD('exec_prob_tab_neigh requires prior particle sampling (in exec_prob_align)')
        endif
        ! PREPARE REFERENCES, SIGMAS, POLAR_CORRCALC, PTCLS
        call prepare_refs_sigmas_ptcls( params, build, cline, tmp_imgs, tmp_imgs_pad, nptcls, params%which_iter,&
                                        do_polar=(params%l_polar .and. (.not.cline%defined('vol1'))) )
        ! Build polar particle images
        call build_batch_particles(params, build, nptcls, pinds, tmp_imgs, tmp_imgs_pad)
        if( .not. allocated(build%subspace_inds) )then
            THROW_HARD('exec_prob_tab_neigh requires neighborhood subspace indices; enable l_neigh and set nspace_sub')
        endif
        if( size(build%subspace_inds) /= params%nspace_sub )then
            THROW_HARD('exec_prob_tab_neigh: size(subspace_inds) must equal nspace_sub')
        endif
        call neigh_map%new(build%subspace_inds, params%nspace_sub)
        call build_neigh_mask_from_subspace_peaks(neigh_mask, ptcl_state)
        call eulprob_obj_part_neigh%new(params, build, pinds, neigh_mask, ptcl_state)
        call eulprob_obj_part_neigh%fill_tab
        call eulprob_obj_part_neigh%ref_assign
        ! Write local assignment map for this part.
        fname = string(ASSIGNMENT_FBODY)//'_neigh_'//int2str_pad(params%part,params%numlen)//'.dat'
        call eulprob_obj_part_neigh%write_assignment(fname)
        call eulprob_obj_part_neigh%kill
        call neigh_map%kill
        if( allocated(neigh_mask) ) deallocate(neigh_mask)
        if( allocated(ptcl_state) ) deallocate(ptcl_state)
        call killimgbatch(build)
        call build%pftc%kill
        call build%kill_general_tbox
        call dealloc_imgarr(tmp_imgs)
        call dealloc_imgarr(tmp_imgs_pad)
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_tab_neigh'))
        call simple_end('**** SIMPLE_PROB_TAB_NEIGH NORMAL STOP ****', print_simple=.false.)

    contains

        subroutine build_neigh_mask_from_subspace_peaks(mask, state_out)
            logical, allocatable, intent(out) :: mask(:,:)
            integer, allocatable, intent(out) :: state_out(:)
            type(pftc_shsrch_grad) :: grad_shsrch_obj(nthr_glob)
            type(ori) :: o_prev
            integer, allocatable :: inpl_sorted_inds(:,:), coarse_rank(:,:), peak_sub_idxs(:,:), full2sub(:)
            real,    allocatable :: inpl_dists(:,:), inpl_dists_sorted(:,:), coarse_best_dist(:,:)
            real,    allocatable :: inpl_athres(:)
            real    :: lims(2,2), lims_init(2,2)
            real    :: cxy(3), cxy_shift(2)
            logical :: do_shift_first
            integer :: nrots, nspace_sub, npeak_use
            integer :: i, isub, ithr, iptcl, istate, iproj_full, irot, iref_start, iref_prev
            logical :: prev_mask(params%nspace)
            nspace_sub     = params%nspace_sub
            npeak_use      = max(1, min(params%npeaks, nspace_sub))
            nrots          = build%pftc%get_nrots()
            do_shift_first = params%l_sh_first .and. params%l_doshift
            allocate(mask(params%nspace, nptcls), source=.false.)
            allocate(state_out(nptcls), source=0)
            allocate(inpl_athres(params%nstates), source=0.)
            do istate = 1, params%nstates
                inpl_athres(istate) = calc_athres(build%spproj_field, 'dist_inpl', params%prob_athres, state=istate)
            enddo
            allocate(inpl_dists(nrots, nthr_glob), inpl_dists_sorted(nrots, nthr_glob), inpl_sorted_inds(nrots, nthr_glob))
            allocate(coarse_best_dist(nspace_sub, nthr_glob), coarse_rank(nspace_sub, nthr_glob), peak_sub_idxs(npeak_use, nthr_glob))
            full2sub = neigh_map%get_full2sub_map()
            if( do_shift_first )then
                lims(:,1)      = -params%trs
                lims(:,2)      =  params%trs
                lims_init(:,1) = -SHC_INPL_TRSHWDTH
                lims_init(:,2) =  SHC_INPL_TRSHWDTH
                do ithr = 1, nthr_glob
                    call grad_shsrch_obj(ithr)%new(build, lims, lims_init=lims_init, shbarrier=params%shbarrier, &
                                                   maxits=params%maxits_sh, opt_angle=.true., coarse_init=.true.)
                enddo
            endif
            !$omp parallel do default(shared) private(i,ithr,iptcl,o_prev,istate,isub,iproj_full,irot,iref_start,iref_prev,cxy,cxy_shift,prev_mask) proc_bind(close) schedule(static)
            do i = 1, nptcls
                iptcl = pinds(i)
                ithr  = omp_get_thread_num() + 1
                call build%spproj_field%get_ori(iptcl, o_prev)
                istate = o_prev%get_state()
                if( istate < 1 .or. istate > params%nstates )then
                    THROW_HARD('exec_prob_tab_neigh: particle state out of range in previous orientation')
                endif
                state_out(i) = istate
                cxy_shift    = [0.,0.]
                iproj_full = build%eulspace%find_closest_proj(o_prev)
                if( do_shift_first .and. iproj_full >= 1 .and. iproj_full <= params%nspace )then
                    iref_start = (istate-1)*params%nspace
                    iref_prev  = iref_start + iproj_full
                    irot       = build%pftc%get_roind(360. - o_prev%e3get())
                    call grad_shsrch_obj(ithr)%set_indices(iref_prev, iptcl)
                    cxy = grad_shsrch_obj(ithr)%minimize(irot=irot, sh_rot=.false.)
                    if( irot /= 0 ) cxy_shift = cxy(2:3)
                endif

                if( istate >= 1 .and. istate <= params%nstates )then
                    do isub = 1, nspace_sub
                        iproj_full = build%subspace_inds(isub)
                        call build%pftc%gen_objfun_vals((istate-1)*params%nspace + iproj_full, iptcl, cxy_shift, inpl_dists(:,ithr))
                        inpl_dists(:,ithr) = eulprob_dist_switch(inpl_dists(:,ithr), params%cc_objfun)
                        irot = angle_sampling(inpl_dists(:,ithr), inpl_dists_sorted(:,ithr), inpl_sorted_inds(:,ithr), inpl_athres(istate), params%prob_athres)
                        coarse_best_dist(isub,ithr) = inpl_dists(irot,ithr)
                        coarse_rank(isub,ithr)      = isub
                    enddo
                    call hpsort(coarse_best_dist(:,ithr), coarse_rank(:,ithr))
                    peak_sub_idxs(:,ithr) = coarse_rank(1:npeak_use,ithr)
                    iproj_full = build%eulspace%find_closest_proj(o_prev)
                    call neigh_map%get_neighbors_mask_pooled(peak_sub_idxs(:,ithr), mask(:,i))
                    if( iproj_full >= 1 .and. iproj_full <= params%nspace )then
                        call neigh_map%get_neighbors_mask(full2sub(iproj_full), prev_mask)
                        mask(:,i) = mask(:,i) .or. prev_mask
                    endif
                else
                    iproj_full = build%eulspace%find_closest_proj(o_prev)
                    if( iproj_full >= 1 .and. iproj_full <= params%nspace )then
                        call neigh_map%get_neighbors_mask(full2sub(iproj_full), mask(:,i))
                    else
                        mask(:,i) = .true.
                    endif
                endif
                call o_prev%kill
            enddo
            !$omp end parallel do
            if( do_shift_first )then
                do ithr = 1, nthr_glob
                    call grad_shsrch_obj(ithr)%kill
                enddo
            endif
            if( allocated(full2sub)          ) deallocate(full2sub)
            if( allocated(inpl_athres)       ) deallocate(inpl_athres)
            if( allocated(inpl_dists)        ) deallocate(inpl_dists)
            if( allocated(inpl_dists_sorted) ) deallocate(inpl_dists_sorted)
            if( allocated(inpl_sorted_inds)  ) deallocate(inpl_sorted_inds)
            if( allocated(coarse_best_dist)  ) deallocate(coarse_best_dist)
            if( allocated(coarse_rank)       ) deallocate(coarse_rank)
            if( allocated(peak_sub_idxs)     ) deallocate(peak_sub_idxs)
        end subroutine build_neigh_mask_from_subspace_peaks

    end subroutine exec_prob_tab_neigh

    subroutine exec_prob_align( self, cline )
        use simple_eul_prob_tab,        only: eul_prob_tab
        use simple_strategy2D3D_common, only: sample_ptcls4fillin, sample_ptcls4update
        use simple_builder,             only: builder
        class(commander_prob_align), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,            allocatable :: pinds(:)
        type(string)                    :: fname
        type(builder)                   :: build
        type(parameters)                :: params
        type(commander_prob_tab)        :: xprob_tab
        type(eul_prob_tab)              :: eulprob_obj_glob
        type(cmdline)                   :: cline_prob_tab
        type(qsys_env)                  :: qenv
        type(chash)                     :: job_descr
        integer :: nptcls, ipart
        call cline%set('mkdir',  'no')
        call cline%set('stream', 'no')
        call build%init_params_and_build_general_tbox(cline, params, do3d=.true.)
        if( params%startit == 1 ) call build%spproj_field%clean_entry('updatecnt', 'sampled')
        ! sampled incremented
        if( params%l_fillin .and. mod(params%startit,5) == 0 )then
            call sample_ptcls4fillin(params, build, [1,params%nptcls], .true., nptcls, pinds)
        else
            call sample_ptcls4update(params, build, [1,params%nptcls], .true., nptcls, pinds)
        endif
        ! communicate to project file
        call build%spproj%write_segment_inside(params%oritype)
        ! more prep
        call eulprob_obj_glob%new(params, build, pinds)
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
        ! reading corrs from all parts
        if( str_has_substr(params%refine, 'prob_state') )then
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
        ! cleanup
        call eulprob_obj_glob%kill
        call cline_prob_tab%kill
        call qenv%kill
        call job_descr%kill
        call build%kill_general_tbox
        call qsys_job_finished(params, string('simple_commanders_refine3D :: exec_prob_align'))
        call qsys_cleanup(params)
        call simple_end('**** SIMPLE_PROB_ALIGN NORMAL STOP ****', print_simple=.false.)
    end subroutine exec_prob_align

end module simple_commanders_prob