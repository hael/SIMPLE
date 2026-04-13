!@descr: neighborhood extension of probabilistic 3D search table.
! The sparse neighborhood is geometrically sparse and searched independently per state.
module simple_eul_prob_tab_neigh
use simple_pftc_srch_api
use simple_builder,            only: builder
use simple_eul_prob_tab,       only: eul_prob_tab, angle_sampling, calc_num2sample, calc_athres, eulprob_dist_switch
use simple_pftc_shsrch_grad,   only: pftc_shsrch_grad
use simple_ori,                only: ori
use simple_eulspace_neigh_map, only: eulspace_neigh_map
use simple_strategy_tree_helpers, only: choose_next_child_prob, INVALID_CORR
implicit none

public :: eul_prob_tab_neigh
private
#include "simple_local_flags.inc"

type, extends(eul_prob_tab) :: eul_prob_tab_neigh
    integer, allocatable    :: eval_touched_refs(:,:)
    integer, allocatable    :: eval_touched_counts(:)
    integer                 :: eval_max_touched = 0
contains
    procedure :: new_neigh
    procedure :: fill_tab      => fill_tab_neigh
    procedure :: ref_normalize => ref_normalize_neigh
    procedure :: ref_assign    => ref_assign_neigh
    procedure :: kill          => kill_neigh
    procedure :: read_tabs_to_glob
end type eul_prob_tab_neigh

contains

    subroutine new_neigh( self, params, build, pinds, empty_okay )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(parameters), target, intent(in)    :: params
        class(builder),    target, intent(in)    :: build
        integer,                   intent(in)    :: pinds(:)
        logical, optional,         intent(in)    :: empty_okay
        integer :: nsubs
        call self%kill
        call self%eul_prob_tab%new(params, build, pinds, empty_okay)
        nsubs = size(self%b_ptr%subspace_inds)
        if( nsubs < 1 )then
            THROW_HARD('simple_eul_prob_tab_neigh::new_neigh; empty subspace indices')
        endif
        ! Sparse support now spans neighborhoods across all active states.
        self%eval_max_touched = max(1, self%nrefs)
        allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        allocate(self%eval_touched_counts(self%nptcls),                     source=0)
    end subroutine new_neigh

    subroutine fill_tab_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        ! Per particle:
        ! 1) get previous-state context and estimate one shift seed
        ! 2) coarse-score one representative per subspace for every active state
        ! 3) pick top subspaces per state and pool neighborhoods per state
        ! 4) evaluate all refs whose projection falls in each state's pooled neighborhood
        ! 5) optionally refine the best evaluated refs with shift search

        type :: coarse_search_ws
            real,    allocatable :: best_subspace_dist(:,:,:)
            real,    allocatable :: best_subspace_dist_work(:,:,:)
            integer, allocatable :: peak_subspace_inds(:,:,:)
            integer, allocatable :: peak_subspace_count(:,:)
            integer, allocatable :: pooled_sub_inds(:,:,:)
            logical, allocatable :: neigh_proj_mask(:,:,:)
        end type coarse_search_ws

        type :: eval_ws
            integer, allocatable :: evaluated_ref_ids(:,:)
            integer, allocatable :: best_eval_locs(:,:)
            integer, allocatable :: fullref_to_sparse_ref(:)
            real,    allocatable :: evaluated_ref_dists(:,:)
            integer, allocatable :: cand_ref_ids(:,:)
            integer, allocatable :: cand_inpl(:,:)
            integer, allocatable :: cand_counts(:)
            real,    allocatable :: cand_dists(:,:)
            real,    allocatable :: cand_x(:,:)
            real,    allocatable :: cand_y(:,:)
            logical, allocatable :: cand_has_sh(:,:)
        end type eval_ws

        type(eulspace_neigh_map) :: neigh_map
        type(pftc_shsrch_grad)   :: grad_shsrch_obj(nthr_glob)
        type(ori)                :: o_prev
        type(coarse_search_ws)   :: coarse_ws
        type(eval_ws)            :: eval_work
        integer, allocatable     :: inds_sorted(:,:)
        real,    allocatable     :: inpl_athres(:), dists_inpl(:,:), dists_inpl_sorted(:,:)
        integer :: i, ri, istate, ithr, max_refs_to_refine
        integer :: iptcl, nsubs, npeak_target, n, si, iref_full
        real    :: lims(2,2), lims_init(2,2), shift_seed(3)
        logical :: l_use_tree_descent
        call seed_rnd
        nsubs = size(self%b_ptr%subspace_inds)
        npeak_target = min(max(1, self%p_ptr%npeaks), nsubs)
        l_use_tree_descent = trim(self%p_ptr%refine) == 'prob_tree'
        call neigh_map%new(self%b_ptr%subspace_full2sub_map, nsubs)
        if( self%eval_max_touched < 1 ) self%eval_max_touched = 1
        if( .not. allocated(self%eval_touched_refs)   ) allocate(self%eval_touched_refs(self%eval_max_touched,self%nptcls), source=0)
        if( .not. allocated(self%eval_touched_counts) ) allocate(self%eval_touched_counts(self%nptcls),                     source=0)
        call clear_sparse_eval_table()
        allocate(inpl_athres(self%p_ptr%nstates), source=self%p_ptr%prob_athres)
        allocate(dists_inpl(self%b_ptr%pftc%get_nrots(),nthr_glob),&
        &dists_inpl_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob),&
        &inds_sorted(self%b_ptr%pftc%get_nrots(),nthr_glob))
        ! determine max number of neighbors to refine per active state
        max_refs_to_refine = 0
        do si = 1, self%nstates
            istate = self%ssinds(si)
            call calc_num2sample(self%b_ptr%spproj_field, self%p_ptr%nspace, 'dist', n, self%p_ptr%prob_athres, state=istate)
            max_refs_to_refine  = max(max_refs_to_refine, n)
            inpl_athres(istate) = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=istate)
        enddo
        if( allocated(eval_work%best_eval_locs) ) deallocate(eval_work%best_eval_locs)
        allocate(eval_work%best_eval_locs(max_refs_to_refine,nthr_glob), eval_work%evaluated_ref_ids(self%nrefs,nthr_glob), source=0)
        allocate(eval_work%fullref_to_sparse_ref(self%p_ptr%nstates*self%p_ptr%nspace),                                     source=0)
        do ri = 1, self%nrefs
            iref_full = (self%sinds(ri)-1)*self%p_ptr%nspace + self%jinds(ri)
            if( iref_full >= 1 .and. iref_full <= size(eval_work%fullref_to_sparse_ref) ) eval_work%fullref_to_sparse_ref(iref_full) = ri
        enddo
        allocate(coarse_ws%neigh_proj_mask(self%p_ptr%nspace,self%p_ptr%nstates,nthr_glob), source=.false.)
        allocate(eval_work%evaluated_ref_dists(self%nrefs,nthr_glob),                   source=huge(1.0))
        allocate(coarse_ws%best_subspace_dist(nsubs,self%p_ptr%nstates,nthr_glob),      source=huge(1.0))
        allocate(coarse_ws%best_subspace_dist_work(nsubs,self%p_ptr%nstates,nthr_glob), source=huge(1.0))
        allocate(coarse_ws%peak_subspace_inds(nsubs,self%p_ptr%nstates,nthr_glob),      source=0)
        allocate(coarse_ws%pooled_sub_inds(npeak_target,self%p_ptr%nstates,nthr_glob),  source=0)
        allocate(coarse_ws%peak_subspace_count(self%p_ptr%nstates,nthr_glob),           source=0)
        allocate(eval_work%cand_ref_ids(self%nrefs,nthr_glob), source=0)
        allocate(eval_work%cand_inpl(self%nrefs,nthr_glob),    source=0)
        allocate(eval_work%cand_counts(nthr_glob),             source=0)
        allocate(eval_work%cand_dists(self%nrefs,nthr_glob),   source=huge(1.0))
        allocate(eval_work%cand_x(self%nrefs,nthr_glob),       source=0.)
        allocate(eval_work%cand_y(self%nrefs,nthr_glob),       source=0.)
        allocate(eval_work%cand_has_sh(self%nrefs,nthr_glob),  source=.false.)
        if( self%p_ptr%l_doshift )then
            ! make shift search objects
            lims(:,1)      = -self%p_ptr%trs
            lims(:,2)      =  self%p_ptr%trs
            lims_init(:,1) = -SHC_INPL_TRSHWDTH
            lims_init(:,2) =  SHC_INPL_TRSHWDTH
            do ithr = 1,nthr_glob
                call grad_shsrch_obj(ithr)%new(self%b_ptr, lims, lims_init=lims_init, shbarrier=self%p_ptr%shbarrier,&
                    &maxits=self%p_ptr%maxits_sh, opt_angle=.true., coarse_init=.true.)
            end do
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                call process_particle(i, iptcl, ithr, .true., shift_seed)
            enddo
            !$omp end parallel do
        else
            ! no shift search - evaluate pooled neighborhoods with zero shift only
            !$omp parallel do default(shared) private(i,iptcl,ithr,shift_seed) proc_bind(close) schedule(static)
            do i = 1, self%nptcls
                iptcl = self%pinds(i)
                ithr  = omp_get_thread_num() + 1
                shift_seed = 0.
                call process_particle(i, iptcl, ithr, .false., shift_seed)
            enddo
            !$omp end parallel do
        endif
        do ithr = 1,nthr_glob
            call grad_shsrch_obj(ithr)%kill
        end do
        call neigh_map%kill
        deallocate(coarse_ws%peak_subspace_count, coarse_ws%pooled_sub_inds, coarse_ws%peak_subspace_inds,&
        &coarse_ws%best_subspace_dist_work, coarse_ws%best_subspace_dist, eval_work%evaluated_ref_dists,&
        &eval_work%evaluated_ref_ids, coarse_ws%neigh_proj_mask, eval_work%best_eval_locs, eval_work%fullref_to_sparse_ref,&
        &eval_work%cand_ref_ids, eval_work%cand_inpl, eval_work%cand_counts, eval_work%cand_dists, eval_work%cand_x,&
        &eval_work%cand_y, eval_work%cand_has_sh,&
        &inds_sorted, dists_inpl_sorted, dists_inpl, inpl_athres)

    contains

        subroutine clear_sparse_eval_table()
            integer :: i_loc, k_loc, ri_loc
            !$omp parallel do default(shared) private(i_loc,k_loc,ri_loc) proc_bind(close) schedule(static)
            do i_loc = 1, self%nptcls
                do k_loc = 1, self%eval_touched_counts(i_loc)
                    ri_loc = self%eval_touched_refs(k_loc,i_loc)
                    if( ri_loc > 0 )then
                        self%loc_tab(ri_loc,i_loc)%dist   = huge(1.0)
                        self%loc_tab(ri_loc,i_loc)%inpl   = 0
                        self%loc_tab(ri_loc,i_loc)%x      = 0.
                        self%loc_tab(ri_loc,i_loc)%y      = 0.
                        self%loc_tab(ri_loc,i_loc)%has_sh = .false.
                    endif
                enddo
                self%eval_touched_counts(i_loc) = 0
            enddo
            !$omp end parallel do
        end subroutine clear_sparse_eval_table

        subroutine process_particle(i_loc, iptcl_loc, ithr_loc, l_with_shift, shift_seed_loc)
            integer, intent(in)    :: i_loc, iptcl_loc, ithr_loc
            logical, intent(in)    :: l_with_shift
            real,    intent(inout) :: shift_seed_loc(3)
            type(ori) :: o_prev_loc
            integer :: prev_state_loc, prev_proj_loc
            integer :: neval_loc
            call get_particle_context(iptcl_loc, o_prev_loc, prev_state_loc, prev_proj_loc)
            call estimate_shift_seed(ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc, o_prev_loc, l_with_shift, shift_seed_loc)
            call find_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
            if( l_use_tree_descent )then
                call evaluate_tree_descent_or_fallback(i_loc, ithr_loc, iptcl_loc, prev_proj_loc, shift_seed_loc, l_with_shift, neval_loc)
            else
                call build_pooled_neighborhood(ithr_loc, prev_proj_loc)
                call evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval_loc)
            endif
            call refine_best_neighbors(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, neval_loc, l_with_shift)
            call o_prev_loc%kill
        end subroutine process_particle

        subroutine evaluate_tree_descent_or_fallback(i_loc, ithr_loc, iptcl_loc, prev_proj_loc, shift_seed_loc, l_with_shift, neval_loc)
            integer, intent(in)  :: i_loc, ithr_loc, iptcl_loc, prev_proj_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            integer, intent(out) :: neval_loc
            call clear_candidate_store(ithr_loc)
            call evaluate_tree_descent(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval_loc)
            if( neval_loc > 0 )then
                call flush_candidate_store(i_loc, ithr_loc, neval_loc)
            else
                call build_pooled_neighborhood(ithr_loc, prev_proj_loc)
                call evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval_loc)
            endif
        end subroutine evaluate_tree_descent_or_fallback

        subroutine evaluate_tree_descent(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval_loc)
            use simple_binary_tree, only: bt_node
            integer, intent(in)  :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in)  :: shift_seed_loc(3)
            logical, intent(in)  :: l_with_shift
            integer, intent(out) :: neval_loc
            type(bt_node) :: node_cur, node_root
            logical, allocatable :: tree_seen(:)
            integer :: si_loc, istate_loc, ntrees_loc, ipeak_loc, isub_loc, iproj_loc, itree_loc
            integer :: inode_loc, inode_next_loc
            real    :: score_root, score_l, score_r
            ntrees_loc = self%b_ptr%block_tree%get_n_trees()
            if( ntrees_loc <= 0 )then
                neval_loc = 0
                return
            endif
            allocate(tree_seen(ntrees_loc), source=.false.)
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                if( .not. self%state_exists(istate_loc) ) cycle
                tree_seen = .false.
                do ipeak_loc = 1, coarse_ws%peak_subspace_count(istate_loc,ithr_loc)
                    isub_loc = coarse_ws%peak_subspace_inds(ipeak_loc,istate_loc,ithr_loc)
                    if( isub_loc < 1 .or. isub_loc > nsubs ) cycle
                    iproj_loc = self%b_ptr%subspace_inds(isub_loc)
                    if( iproj_loc < 1 .or. iproj_loc > self%p_ptr%nspace ) cycle
                    itree_loc = self%b_ptr%subspace_full2sub_map(iproj_loc)
                    if( itree_loc < 1 .or. itree_loc > ntrees_loc ) cycle
                    if( tree_seen(itree_loc) ) cycle
                    tree_seen(itree_loc) = .true.
                    node_root = self%b_ptr%block_tree%get_root_node(itree_loc)
                    call eval_tree_ref_fixed_state_local(node_root%ref_idx, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_root)
                    inode_loc = node_root%node_idx
                    do
                        if( inode_loc == 0 ) exit
                        if( self%b_ptr%block_tree%is_leaf(itree_loc, inode_loc) ) exit
                        node_cur = self%b_ptr%block_tree%get_node(itree_loc, inode_loc)
                        if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
                        call eval_child_fixed_state_local(itree_loc, node_cur%left_idx, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_l)
                        call eval_child_fixed_state_local(itree_loc, node_cur%right_idx, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_r)
                        inode_next_loc = choose_next_child_prob(node_cur%left_idx, node_cur%right_idx, score_l, score_r)
                        if( inode_next_loc == 0 ) exit
                        inode_loc = inode_next_loc
                    enddo
                enddo
            enddo
            neval_loc = eval_work%cand_counts(ithr_loc)
            deallocate(tree_seen)
        end subroutine evaluate_tree_descent

        subroutine eval_child_fixed_state_local(itree_loc, child_idx_loc, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_loc)
            use simple_binary_tree, only: bt_node
            integer, intent(in) :: itree_loc, child_idx_loc, istate_loc, i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            real,    intent(out) :: score_loc
            type(bt_node) :: node_child
            score_loc = INVALID_CORR
            if( child_idx_loc == 0 ) return
            node_child = self%b_ptr%block_tree%get_node(itree_loc, child_idx_loc)
            if( node_child%ref_idx == 0 ) return
            call eval_tree_ref_fixed_state_local(node_child%ref_idx, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_loc)
        end subroutine eval_child_fixed_state_local

        subroutine eval_tree_ref_fixed_state_local(ref_idx_loc, istate_loc, i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, score_loc)
            integer, intent(in) :: ref_idx_loc, istate_loc, i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            real,    intent(out) :: score_loc
            integer :: iref_full_loc, irot_loc, ri_loc
            real    :: rotmat_loc(2,2), rotated_shift_loc(2), dist_best
            score_loc = INVALID_CORR
            if( ref_idx_loc < 1 .or. ref_idx_loc > self%p_ptr%nspace ) return
            if( .not. self%state_exists(istate_loc) ) return
            if( .not. self%proj_exists(ref_idx_loc,istate_loc) ) return
            iref_full_loc = (istate_loc-1)*self%p_ptr%nspace + ref_idx_loc
            if( l_with_shift )then
                call self%b_ptr%pftc%gen_objfun_vals(iref_full_loc, iptcl_loc, shift_seed_loc(2:3), dists_inpl(:,ithr_loc))
            else
                call self%b_ptr%pftc%gen_objfun_vals(iref_full_loc, iptcl_loc, [0.,0.], dists_inpl(:,ithr_loc))
            endif
            dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
            irot_loc = angle_sampling(dists_inpl(:,ithr_loc), dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc), inpl_athres(istate_loc), self%p_ptr%prob_athres)
            dist_best = dists_inpl(irot_loc,ithr_loc)
            ri_loc = eval_work%fullref_to_sparse_ref(iref_full_loc)
            if( ri_loc > 0 )then
                if( l_with_shift )then
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                    rotated_shift_loc = matmul(shift_seed_loc(2:3), rotmat_loc)
                    call store_solution_local(ithr_loc, ri_loc, dist_best, irot_loc, rotated_shift_loc(1), rotated_shift_loc(2), .true.)
                else
                    call store_solution_local(ithr_loc, ri_loc, dist_best, irot_loc, 0., 0., .false.)
                endif
            endif
            score_loc = -dist_best
        end subroutine eval_tree_ref_fixed_state_local

        subroutine get_particle_context(iptcl_loc, o_prev_loc, prev_state_loc, prev_proj_loc)
            integer,   intent(in)    :: iptcl_loc
            type(ori), intent(inout) :: o_prev_loc
            integer,   intent(out)   :: prev_state_loc, prev_proj_loc
            call self%b_ptr%spproj_field%get_ori(iptcl_loc, o_prev_loc)
            prev_state_loc = o_prev_loc%get_state()
            prev_proj_loc  = self%b_ptr%eulspace%find_closest_proj(o_prev_loc)
        end subroutine get_particle_context

        subroutine estimate_shift_seed(ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc, o_prev_loc, l_with_shift, shift_seed_loc)
            integer,   intent(in)    :: ithr_loc, iptcl_loc, prev_state_loc, prev_proj_loc
            type(ori), intent(inout) :: o_prev_loc
            logical,   intent(in)    :: l_with_shift
            real,      intent(inout) :: shift_seed_loc(3)
            integer :: irot_loc, iref_start_loc
            if( .not. l_with_shift )then
                shift_seed_loc = 0.
                return
            endif
            irot_loc = self%b_ptr%pftc%get_roind(360.-o_prev_loc%e3get())
            if( self%state_exists(prev_state_loc) .and. self%proj_exists(prev_proj_loc,prev_state_loc) )then
                iref_start_loc = (prev_state_loc-1)*self%p_ptr%nspace
                call grad_shsrch_obj(ithr_loc)%set_indices(iref_start_loc + prev_proj_loc, iptcl_loc)
                shift_seed_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.false.)
                if( irot_loc == 0 ) shift_seed_loc(2:3) = 0.
            else
                shift_seed_loc = 0.
            endif
        end subroutine estimate_shift_seed

        subroutine find_peak_subspaces(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            integer :: si_loc, istate_loc, isub_loc, full_ref_subspace_loc, irot_loc, ri_loc, ipeak_loc, coarse_proj_loc
            real    :: rotmat_loc(2,2), rotated_shift_loc(2)
            coarse_ws%best_subspace_dist(:,:,ithr_loc) = huge(1.0)
            coarse_ws%peak_subspace_count(:,ithr_loc)  = 0
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                if( .not. self%state_exists(istate_loc) ) cycle
                do isub_loc = 1, nsubs
                    coarse_proj_loc = self%b_ptr%subspace_inds(isub_loc)
                    if( .not. self%proj_exists(coarse_proj_loc, istate_loc) ) cycle
                    full_ref_subspace_loc = (istate_loc-1)*self%p_ptr%nspace + coarse_proj_loc
                    if( l_with_shift )then
                        call self%b_ptr%pftc%gen_objfun_vals(full_ref_subspace_loc, iptcl_loc, shift_seed_loc(2:3), dists_inpl(:,ithr_loc))
                    else
                        call self%b_ptr%pftc%gen_objfun_vals(full_ref_subspace_loc, iptcl_loc, [0.,0.], dists_inpl(:,ithr_loc))
                    endif
                    dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                    irot_loc = minloc(dists_inpl(:,ithr_loc), dim=1)
                    coarse_ws%best_subspace_dist(isub_loc,istate_loc,ithr_loc) = dists_inpl(irot_loc,ithr_loc)
                    ri_loc = eval_work%fullref_to_sparse_ref(full_ref_subspace_loc)
                    if( ri_loc > 0 )then
                        if( l_with_shift )then
                            call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                            rotated_shift_loc = matmul(shift_seed_loc(2:3), rotmat_loc)
                            call record_sparse_eval(i_loc, ri_loc, dists_inpl(irot_loc,ithr_loc), irot_loc, rotated_shift_loc(1), rotated_shift_loc(2), .true.)
                        else
                            call record_sparse_eval(i_loc, ri_loc, dists_inpl(irot_loc,ithr_loc), irot_loc, 0., 0., .false.)
                        endif
                    endif
                enddo
                coarse_ws%best_subspace_dist_work(:,istate_loc,ithr_loc) = coarse_ws%best_subspace_dist(:,istate_loc,ithr_loc)
                do ipeak_loc = 1, npeak_target
                    coarse_ws%peak_subspace_inds(ipeak_loc,istate_loc,ithr_loc) = minloc(coarse_ws%best_subspace_dist_work(:,istate_loc,ithr_loc), dim=1)
                    isub_loc = coarse_ws%peak_subspace_inds(ipeak_loc,istate_loc,ithr_loc)
                    if( coarse_ws%best_subspace_dist_work(isub_loc,istate_loc,ithr_loc) >= huge(1.0) )then
                        coarse_ws%peak_subspace_inds(ipeak_loc,istate_loc,ithr_loc) = 0
                        exit
                    endif
                    coarse_ws%best_subspace_dist_work(isub_loc,istate_loc,ithr_loc) = huge(1.0)
                    coarse_ws%peak_subspace_count(istate_loc,ithr_loc) = coarse_ws%peak_subspace_count(istate_loc,ithr_loc) + 1
                enddo
            enddo
        end subroutine find_peak_subspaces

        subroutine build_pooled_neighborhood(ithr_loc, prev_proj_loc)
            integer, intent(in) :: ithr_loc, prev_proj_loc
            integer :: si_loc, istate_loc, npeak_found_loc, isub_loc, coarse_proj_loc, iproj_loc
            coarse_ws%neigh_proj_mask(:,:,ithr_loc) = .false.
            do si_loc = 1, self%nstates
                istate_loc = self%ssinds(si_loc)
                npeak_found_loc = coarse_ws%peak_subspace_count(istate_loc,ithr_loc)
                if( npeak_found_loc > 0 )then
                    coarse_ws%pooled_sub_inds(1:npeak_found_loc,istate_loc,ithr_loc) = coarse_ws%peak_subspace_inds(1:npeak_found_loc,istate_loc,ithr_loc)
                    call neigh_map%get_neighbors_mask_pooled(coarse_ws%pooled_sub_inds(1:npeak_found_loc,istate_loc,ithr_loc), coarse_ws%neigh_proj_mask(:,istate_loc,ithr_loc))
                else
                    coarse_proj_loc = max(1, min(self%p_ptr%nspace, prev_proj_loc))
                    if( .not. self%proj_exists(coarse_proj_loc,istate_loc) )then
                        coarse_proj_loc = 0
                        do iproj_loc = 1, self%p_ptr%nspace
                            if( self%proj_exists(iproj_loc,istate_loc) )then
                                coarse_proj_loc = iproj_loc
                                exit
                            endif
                        enddo
                    endif
                    if( coarse_proj_loc < 1 ) cycle
                    isub_loc = self%b_ptr%subspace_full2sub_map(coarse_proj_loc)
                    if( isub_loc < 1 .or. isub_loc > nsubs ) isub_loc = 1
                    call neigh_map%get_neighbors_mask(isub_loc, coarse_ws%neigh_proj_mask(:,istate_loc,ithr_loc))
                endif
            enddo
        end subroutine build_pooled_neighborhood

        subroutine evaluate_neighborhood(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, l_with_shift, neval_loc)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            integer, intent(out) :: neval_loc
            integer :: ri_loc, istate_loc, iproj_loc, irot_loc, iref_loc
            real    :: rotmat_loc(2,2), rotated_shift_loc(2)
            call clear_candidate_store(ithr_loc)
            do ri_loc = 1, self%nrefs
                istate_loc = self%sinds(ri_loc)
                iproj_loc  = self%jinds(ri_loc)
                if( .not. coarse_ws%neigh_proj_mask(iproj_loc,istate_loc,ithr_loc) ) cycle
                iref_loc = (istate_loc-1)*self%p_ptr%nspace + iproj_loc
                if( l_with_shift )then
                    call self%b_ptr%pftc%gen_objfun_vals(iref_loc, iptcl_loc, shift_seed_loc(2:3), dists_inpl(:,ithr_loc))
                else
                    call self%b_ptr%pftc%gen_objfun_vals(iref_loc, iptcl_loc, [0.,0.], dists_inpl(:,ithr_loc))
                endif
                dists_inpl(:,ithr_loc) = eulprob_dist_switch(dists_inpl(:,ithr_loc), self%p_ptr%cc_objfun)
                irot_loc = angle_sampling(dists_inpl(:,ithr_loc), dists_inpl_sorted(:,ithr_loc), inds_sorted(:,ithr_loc), inpl_athres(istate_loc), self%p_ptr%prob_athres)
                if( l_with_shift )then
                    call rotmat2d(self%b_ptr%pftc%get_rot(irot_loc), rotmat_loc)
                    rotated_shift_loc = matmul(shift_seed_loc(2:3), rotmat_loc)
                    call store_solution_local(ithr_loc, ri_loc, dists_inpl(irot_loc,ithr_loc), irot_loc, rotated_shift_loc(1), rotated_shift_loc(2), .true.)
                else
                    call store_solution_local(ithr_loc, ri_loc, dists_inpl(irot_loc,ithr_loc), irot_loc, 0., 0., .false.)
                endif
            enddo
            call flush_candidate_store(i_loc, ithr_loc, neval_loc)
        end subroutine evaluate_neighborhood

        subroutine clear_candidate_store(ithr_loc)
            integer, intent(in) :: ithr_loc
            eval_work%cand_counts(ithr_loc) = 0
        end subroutine clear_candidate_store

        subroutine store_solution_local(ithr_loc, ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc)
            integer, intent(in) :: ithr_loc, ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            integer :: j_loc, nloc
            nloc = eval_work%cand_counts(ithr_loc)
            do j_loc = 1, nloc
                if( eval_work%cand_ref_ids(j_loc,ithr_loc) /= ri_loc ) cycle
                if( dist_loc >= eval_work%cand_dists(j_loc,ithr_loc) ) return
                eval_work%cand_dists(j_loc,ithr_loc)  = dist_loc
                eval_work%cand_inpl(j_loc,ithr_loc)   = irot_loc
                eval_work%cand_x(j_loc,ithr_loc)      = x_loc
                eval_work%cand_y(j_loc,ithr_loc)      = y_loc
                eval_work%cand_has_sh(j_loc,ithr_loc) = has_sh_loc
                return
            enddo
            nloc = nloc + 1
            if( nloc > self%nrefs )then
                THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh; local candidate overflow')
            endif
            eval_work%cand_counts(ithr_loc)         = nloc
            eval_work%cand_ref_ids(nloc,ithr_loc)   = ri_loc
            eval_work%cand_dists(nloc,ithr_loc)     = dist_loc
            eval_work%cand_inpl(nloc,ithr_loc)      = irot_loc
            eval_work%cand_x(nloc,ithr_loc)         = x_loc
            eval_work%cand_y(nloc,ithr_loc)         = y_loc
            eval_work%cand_has_sh(nloc,ithr_loc)    = has_sh_loc
        end subroutine store_solution_local

        subroutine flush_candidate_store(i_loc, ithr_loc, neval_loc)
            integer, intent(in)  :: i_loc, ithr_loc
            integer, intent(out) :: neval_loc
            integer :: j_loc, ri_loc
            neval_loc = eval_work%cand_counts(ithr_loc)
            do j_loc = 1, neval_loc
                ri_loc = eval_work%cand_ref_ids(j_loc,ithr_loc)
                call record_sparse_eval(i_loc, ri_loc, eval_work%cand_dists(j_loc,ithr_loc), eval_work%cand_inpl(j_loc,ithr_loc),&
                &eval_work%cand_x(j_loc,ithr_loc), eval_work%cand_y(j_loc,ithr_loc), eval_work%cand_has_sh(j_loc,ithr_loc))
                eval_work%evaluated_ref_ids(j_loc,ithr_loc)   = ri_loc
                eval_work%evaluated_ref_dists(j_loc,ithr_loc) = eval_work%cand_dists(j_loc,ithr_loc)
            enddo
            if( neval_loc < self%nrefs )then
                eval_work%evaluated_ref_ids(neval_loc+1:self%nrefs,ithr_loc)   = 0
                eval_work%evaluated_ref_dists(neval_loc+1:self%nrefs,ithr_loc) = huge(1.0)
            endif
        end subroutine flush_candidate_store

        subroutine refine_best_neighbors(i_loc, ithr_loc, iptcl_loc, shift_seed_loc, neval_loc, l_with_shift)
            integer, intent(in) :: i_loc, ithr_loc, iptcl_loc, neval_loc
            real,    intent(in) :: shift_seed_loc(3)
            logical, intent(in) :: l_with_shift
            integer :: nrefs_to_refine_loc, j_loc, eval_slot_loc, ri_loc, istate_loc, iproj_loc, irot_loc
            real    :: refined_shift_loc(3)
            if( .not. l_with_shift ) return
            eval_work%best_eval_locs(:,ithr_loc) = 0
            if( neval_loc > 0 )then
                nrefs_to_refine_loc = min(max_refs_to_refine, neval_loc)
                eval_work%best_eval_locs(1:nrefs_to_refine_loc,ithr_loc) = minnloc(eval_work%evaluated_ref_dists(1:neval_loc,ithr_loc), nrefs_to_refine_loc)
            else
                nrefs_to_refine_loc = 0
            endif
            do j_loc = 1, nrefs_to_refine_loc
                eval_slot_loc = eval_work%best_eval_locs(j_loc,ithr_loc)
                if( eval_slot_loc < 1 ) cycle
                ri_loc = eval_work%evaluated_ref_ids(eval_slot_loc,ithr_loc)
                if( ri_loc < 1 ) cycle
                istate_loc = self%sinds(ri_loc)
                iproj_loc  = self%jinds(ri_loc)
                call grad_shsrch_obj(ithr_loc)%set_indices((istate_loc-1)*self%p_ptr%nspace + iproj_loc, iptcl_loc)
                irot_loc = self%loc_tab(ri_loc,i_loc)%inpl
                refined_shift_loc = grad_shsrch_obj(ithr_loc)%minimize(irot=irot_loc, sh_rot=.true., xy_in=shift_seed_loc(2:3))
                if( irot_loc > 0 )then
                    call record_sparse_eval(i_loc, ri_loc, eulprob_dist_switch(refined_shift_loc(1), self%p_ptr%cc_objfun), irot_loc, refined_shift_loc(2), refined_shift_loc(3), .true.)
                endif
            enddo
        end subroutine refine_best_neighbors

        subroutine record_sparse_eval(iptcl_loc, ri_loc, dist_loc, irot_loc, x_loc, y_loc, has_sh_loc)
            integer, intent(in) :: iptcl_loc, ri_loc, irot_loc
            real,    intent(in) :: dist_loc, x_loc, y_loc
            logical, intent(in) :: has_sh_loc
            self%loc_tab(ri_loc,iptcl_loc)%dist   = dist_loc
            self%loc_tab(ri_loc,iptcl_loc)%inpl   = irot_loc
            self%loc_tab(ri_loc,iptcl_loc)%x      = x_loc
            self%loc_tab(ri_loc,iptcl_loc)%y      = y_loc
            self%loc_tab(ri_loc,iptcl_loc)%has_sh = has_sh_loc
            call mark_ref_touched(iptcl_loc, ri_loc)
        end subroutine record_sparse_eval

        subroutine mark_ref_touched(iptcl_loc, ri_loc)
            integer, intent(in) :: iptcl_loc, ri_loc
            integer :: kt
            do kt = 1, self%eval_touched_counts(iptcl_loc)
                if( self%eval_touched_refs(kt,iptcl_loc) == ri_loc ) return
            enddo
            if( self%eval_touched_counts(iptcl_loc) >= self%eval_max_touched )then
                THROW_HARD('simple_eul_prob_tab_neigh::fill_tab_neigh; eval_touched overflow')
            endif
            self%eval_touched_counts(iptcl_loc) = self%eval_touched_counts(iptcl_loc) + 1
            self%eval_touched_refs(self%eval_touched_counts(iptcl_loc),iptcl_loc) = ri_loc
        end subroutine mark_ref_touched

    end subroutine fill_tab_neigh

    subroutine read_tabs_to_glob( self, fbody, nparts, numlen )
        class(eul_prob_tab_neigh), intent(inout) :: self
        class(string),             intent(in)    :: fbody
        integer,                   intent(in)    :: nparts, numlen
        type(string) :: fname
        integer, allocatable :: touched_counts(:)
        integer :: ipart, i, ri, max_touched_loaded, cap
        do ipart = 1, nparts
            fname = fbody//int2str_pad(ipart,numlen)//'.dat'
            call self%read_tab_to_glob(fname)
        enddo
        ! Rebuild sparse bookkeeping from loaded loc_tab so assignment works after write/read cycles.
        allocate(touched_counts(self%nptcls), source=0)
        !$omp parallel do default(shared) private(i,ri) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            do ri = 1, self%nrefs
                if( self%loc_tab(ri,i)%inpl > 0 )then
                    touched_counts(i) = touched_counts(i) + 1
                endif
            enddo
        enddo
        !$omp end parallel do
        max_touched_loaded = max(1, maxval(touched_counts))
        if( .not. allocated(self%eval_touched_refs) )then
            allocate(self%eval_touched_refs(max_touched_loaded,self%nptcls), source=0)
        else
            cap = size(self%eval_touched_refs,1)
            if( cap < max_touched_loaded .or. size(self%eval_touched_refs,2) /= self%nptcls )then
                deallocate(self%eval_touched_refs)
                allocate(self%eval_touched_refs(max_touched_loaded,self%nptcls), source=0)
            else
                self%eval_touched_refs = 0
            endif
        endif
        if( .not. allocated(self%eval_touched_counts) )then
            allocate(self%eval_touched_counts(self%nptcls), source=0)
        else if( size(self%eval_touched_counts) /= self%nptcls )then
            deallocate(self%eval_touched_counts)
            allocate(self%eval_touched_counts(self%nptcls), source=0)
        else
            self%eval_touched_counts = 0
        endif
        self%eval_max_touched = size(self%eval_touched_refs,1)
        !$omp parallel do default(shared) private(i,ri) proc_bind(close) schedule(static)
        do i = 1, self%nptcls
            do ri = 1, self%nrefs
                if( self%loc_tab(ri,i)%inpl > 0 )then
                    self%eval_touched_counts(i) = self%eval_touched_counts(i) + 1
                    self%eval_touched_refs(self%eval_touched_counts(i),i) = ri
                endif
            enddo
        enddo
        !$omp end parallel do
        deallocate(touched_counts)
    end subroutine read_tabs_to_glob

    subroutine ref_normalize_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        real    :: sum_dist_all, min_dist, max_dist
        integer :: i, iref, neval
        ! normalize only over evaluated refs (inpl > 0)
        !$omp parallel do default(shared) proc_bind(close) schedule(static) private(i,sum_dist_all,iref,neval)
        do i = 1, self%nptcls
            sum_dist_all = sum(self%loc_tab(:,i)%dist, mask=(self%loc_tab(:,i)%inpl > 0))
            if( sum_dist_all < TINY )then
                ! dense-style behavior: collapse to zero; tie handling randomizes later.
                neval = count(self%loc_tab(:,i)%inpl > 0)
                if( neval > 0 )then
                    do iref = 1, self%nrefs
                        if( self%loc_tab(iref,i)%inpl > 0 ) self%loc_tab(iref,i)%dist = 0.
                    enddo
                endif
            else
                ! divide only evaluated refs
                do iref = 1, self%nrefs
                    if( self%loc_tab(iref,i)%inpl > 0 )&
                    &self%loc_tab(iref,i)%dist = self%loc_tab(iref,i)%dist / sum_dist_all
                enddo
            endif
        enddo
        !$omp end parallel do
        if( .not. any(self%loc_tab(:,:)%inpl > 0) ) return
        ! min/max normalization over evaluated refs only
        min_dist = minval(self%loc_tab(:,:)%dist, mask=(self%loc_tab(:,:)%inpl > 0))
        max_dist = maxval(self%loc_tab(:,:)%dist, mask=(self%loc_tab(:,:)%inpl > 0))
        if( (max_dist - min_dist) < TINY )then
            THROW_WARN('WARNING: numerical unstability in eul_prob_tab_neigh normalize')
            ! dense-style stochastic tie break over evaluated entries only.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(iref,i)
            do iref = 1, self%nrefs
                do i = 1, self%nptcls
                    if( self%loc_tab(iref,i)%inpl > 0 ) self%loc_tab(iref,i)%dist = ran3()
                enddo
            enddo
            !$omp end parallel do
        else
            where( self%loc_tab(:,:)%inpl > 0 )
                self%loc_tab(:,:)%dist = (self%loc_tab(:,:)%dist - min_dist) / (max_dist - min_dist)
            endwhere
        endif
    end subroutine ref_normalize_neigh

    ! sparse graph traversal on-the-fly for reference assignment
    ! with robust fallback to previous orientation if sparse graph leaves particles without valid candidates
    subroutine ref_assign_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        ! Sparse global assignment:
        ! 1) ensure every particle has at least one evaluated candidate
        ! 2) normalize sparse scores
        ! 3) build ref->particle and particle->ref adjacency
        ! 4) repeatedly assign the best currently available particle per active ref
        ! 5) fall back to best evaluated sparse ref for leftovers

        type :: assign_graph_ws
            integer, allocatable :: ref_counts(:), ref_offsets(:), ref_fill(:), ref_list(:), ref_pos(:), active_refs(:)
            integer, allocatable :: ptcl_counts(:), ptcl_offsets(:), ptcl_fill(:), ptcl_refs(:)
            real,    allocatable :: ref_dists(:)
        end type assign_graph_ws

        type :: assign_frontier_ws
            integer, allocatable :: inds_sorted(:), order(:), work_ptcl(:), sel_refs(:), sel_pos(:)
            real,    allocatable :: iref_dist(:), dists_sorted(:), work_d(:), sel_dists(:), sel_dists_sorted(:)
            logical, allocatable :: ptcl_avail(:)
        end type assign_frontier_ws

        type(assign_graph_ws)    :: graph
        type(assign_frontier_ws) :: frontier
        real,    allocatable :: dists_inpl(:), dists_inpl_sorted(:)
        integer, allocatable :: inds_sorted(:)
        integer   :: i, iref, ri, assigned_iref, assigned_ptcl, istate, fallback_ref
        integer   :: k, idx, nactive, total, m, start, maxref, nleft, assigned_idx, nsel, pos, last_ref
        real      :: projs_athres
        real      :: huge_val
        huge_val = huge(1.0)
        allocate(dists_inpl(self%b_ptr%pftc%get_nrots()), dists_inpl_sorted(self%b_ptr%pftc%get_nrots()), inds_sorted(self%b_ptr%pftc%get_nrots()))
        call seed_empty_ptcls_from_prev_assign()
        call self%ref_normalize()
        call build_sparse_assignment_graph()
        call assign_particles_globally()
        call assign_remaining_particles_from_best_touched_ref()
        deallocate(inds_sorted, dists_inpl_sorted, dists_inpl)

    contains

        subroutine seed_empty_ptcls_from_prev_assign()
            do i = 1, self%nptcls
                call seed_fallback_if_empty(i)
            enddo
        end subroutine seed_empty_ptcls_from_prev_assign

        subroutine build_sparse_assignment_graph()
            allocate(graph%ref_counts(self%nrefs), graph%ref_offsets(self%nrefs+1), graph%ref_fill(self%nrefs), graph%ref_pos(self%nrefs),&
            &frontier%iref_dist(self%nrefs), frontier%dists_sorted(self%nrefs), frontier%inds_sorted(self%nrefs), frontier%ptcl_avail(self%nptcls),&
            &graph%ptcl_counts(self%nptcls), graph%ptcl_offsets(self%nptcls+1), graph%ptcl_fill(self%nptcls))
            graph%ref_counts = 0
            graph%ptcl_counts = 0
            do i = 1, self%nptcls
                do k = 1, self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs       ) cycle
                    if( self%loc_tab(iref,i)%inpl <= 0        ) cycle
                    graph%ref_counts(iref) = graph%ref_counts(iref) + 1
                    graph%ptcl_counts(i)   = graph%ptcl_counts(i) + 1
                enddo
            enddo
            nactive = count(graph%ref_counts > 0)
            allocate(graph%active_refs(max(1,nactive)), source=0)
            if( nactive > 0 ) graph%active_refs(1:nactive) = pack((/(iref, iref=1,self%nrefs)/), graph%ref_counts > 0)
            graph%ref_offsets(1) = 1
            do iref = 1, self%nrefs
                graph%ref_offsets(iref+1) = graph%ref_offsets(iref) + graph%ref_counts(iref)
            enddo
            total = graph%ref_offsets(self%nrefs+1) - 1
            allocate(graph%ref_list(max(1,total)), graph%ref_dists(max(1,total)))
            graph%ptcl_offsets(1) = 1
            do i = 1, self%nptcls
                graph%ptcl_offsets(i+1) = graph%ptcl_offsets(i) + graph%ptcl_counts(i)
            enddo
            allocate(graph%ptcl_refs(max(1,graph%ptcl_offsets(self%nptcls+1)-1)))
            graph%ref_fill = graph%ref_offsets(1:self%nrefs)
            graph%ptcl_fill = graph%ptcl_offsets(1:self%nptcls)
            do i = 1, self%nptcls
                do k = 1, self%eval_touched_counts(i)
                    iref = self%eval_touched_refs(k,i)
                    if( iref < 1 .or. iref > self%nrefs ) cycle
                    if( self%loc_tab(iref,i)%inpl <= 0   ) cycle
                    idx = graph%ref_fill(iref)
                    graph%ref_list(idx)  = i
                    graph%ref_dists(idx) = self%loc_tab(iref,i)%dist
                    graph%ref_fill(iref) = idx + 1
                    idx = graph%ptcl_fill(i)
                    graph%ptcl_refs(idx) = iref
                    graph%ptcl_fill(i)   = idx + 1
                enddo
            enddo
            maxref = max(1, maxval(graph%ref_counts))
            allocate(frontier%work_d(maxref), frontier%order(maxref), frontier%work_ptcl(maxref))
            allocate(frontier%sel_refs(max(1,nactive)), frontier%sel_pos(self%nrefs), frontier%sel_dists(max(1,nactive)),&
            &frontier%sel_dists_sorted(max(1,nactive)))
            frontier%sel_pos = 0
            do idx = 1, nactive
                iref = graph%active_refs(idx)
                m    = graph%ref_counts(iref)
                if( m <= 1 ) cycle
                start = graph%ref_offsets(iref)
                frontier%work_d(1:m)    = graph%ref_dists(start:start+m-1)
                frontier%work_ptcl(1:m) = graph%ref_list(start:start+m-1)
                frontier%order(1:m)     = (/(k,k=1,m)/)
                call hpsort(frontier%work_d(1:m), frontier%order(1:m))
                do k = 1, m
                    graph%ref_dists(start+k-1) = frontier%work_d(k)
                    graph%ref_list(start+k-1)  = frontier%work_ptcl(frontier%order(k))
                enddo
            enddo
        end subroutine build_sparse_assignment_graph

        subroutine assign_particles_globally()
            projs_athres = 0.
            do istate = 1, self%nstates
                projs_athres = max(projs_athres, calc_athres(self%b_ptr%spproj_field, 'dist', self%p_ptr%prob_athres, state=istate))
            enddo
            graph%ref_pos       = 1
            frontier%iref_dist  = huge_val
            frontier%ptcl_avail = .true.
            nleft = self%nptcls
            nsel  = 0
            do idx = 1, nactive
                iref = graph%active_refs(idx)
                call advance_ref_head(iref)
                call sync_frontier_ref(iref)
            enddo
            do while( nleft > 0 )
                if( nsel == 0 ) exit
                assigned_idx = angle_sampling(frontier%sel_dists(1:nsel), frontier%sel_dists_sorted(1:nsel), frontier%inds_sorted(1:nsel), projs_athres, self%p_ptr%prob_athres)
                assigned_iref = frontier%sel_refs(assigned_idx)
                assigned_ptcl = graph%ref_list(graph%ref_offsets(assigned_iref) + graph%ref_pos(assigned_iref) - 1)
                frontier%ptcl_avail(assigned_ptcl) = .false.
                nleft = nleft - 1
                self%assgn_map(assigned_ptcl) = self%loc_tab(assigned_iref,assigned_ptcl)
                do idx = graph%ptcl_offsets(assigned_ptcl), graph%ptcl_offsets(assigned_ptcl+1)-1
                    iref = graph%ptcl_refs(idx)
                    m    = graph%ref_counts(iref)
                    if( graph%ref_pos(iref) > m ) cycle
                    start = graph%ref_offsets(iref)
                    if( graph%ref_list(start + graph%ref_pos(iref) - 1) /= assigned_ptcl ) cycle
                    call advance_ref_head(iref)
                    call sync_frontier_ref(iref)
                enddo
            enddo
        end subroutine assign_particles_globally

        subroutine assign_remaining_particles_from_best_touched_ref()
            do i = 1, self%nptcls
                if( .not. frontier%ptcl_avail(i) ) cycle
                fallback_ref = pick_best_evaluated_ref(i)
                if( fallback_ref == 0 )then
                    call seed_fallback_if_empty(i)
                    fallback_ref = pick_best_evaluated_ref(i)
                endif
                if( fallback_ref == 0 ) fallback_ref = 1
                self%assgn_map(i) = self%loc_tab(fallback_ref,i)
            enddo
        end subroutine assign_remaining_particles_from_best_touched_ref

        integer function pick_best_evaluated_ref(iptcl_loc) result(ri_best)
            integer, intent(in) :: iptcl_loc
            integer :: kt, ri_loc
            real    :: best_dist
            ri_best = 0
            best_dist = huge(1.0)
            do kt = 1, self%eval_touched_counts(iptcl_loc)
                ri_loc = self%eval_touched_refs(kt,iptcl_loc)
                if( ri_loc < 1 .or. ri_loc > self%nrefs ) cycle
                if( self%loc_tab(ri_loc,iptcl_loc)%inpl <= 0 ) cycle
                if( self%loc_tab(ri_loc,iptcl_loc)%dist < best_dist )then
                    best_dist = self%loc_tab(ri_loc,iptcl_loc)%dist
                    ri_best   = ri_loc
                endif
            enddo
        end function pick_best_evaluated_ref

        subroutine mark_ref_touched(iptcl_loc, ri_loc)
            integer, intent(in) :: iptcl_loc, ri_loc
            integer :: kt
            do kt = 1, self%eval_touched_counts(iptcl_loc)
                if( self%eval_touched_refs(kt,iptcl_loc) == ri_loc ) return
            enddo
            if( self%eval_touched_counts(iptcl_loc) >= self%eval_max_touched )then
                THROW_HARD('simple_eul_prob_tab_neigh::ref_assign_neigh; eval_touched overflow')
            endif
            self%eval_touched_counts(iptcl_loc) = self%eval_touched_counts(iptcl_loc) + 1
            self%eval_touched_refs(self%eval_touched_counts(iptcl_loc),iptcl_loc) = ri_loc
        end subroutine mark_ref_touched

        subroutine seed_fallback_if_empty(iptcl_loc)
            integer, intent(in) :: iptcl_loc
            type(ori) :: o_prev_loc
            integer :: istate_loc, iproj_loc, irot_loc, fallback_state, fallback_proj, fallback_ref_full
            integer :: ri_loc, fallback_ref_loc
            real    :: inpl_athres_state, sh_seed(2)
            if( pick_best_evaluated_ref(iptcl_loc) > 0 ) return
            call self%b_ptr%spproj_field%get_ori(self%pinds(iptcl_loc), o_prev_loc)
            istate_loc = o_prev_loc%get_state()
            if( istate_loc < 1 .or. istate_loc > self%p_ptr%nstates ) istate_loc = 1
            iproj_loc = self%b_ptr%eulspace%find_closest_proj(o_prev_loc)
            iproj_loc = max(1, min(self%p_ptr%nspace, iproj_loc))
            fallback_ref_loc = 0
            do ri_loc = 1, self%nrefs
                if( self%sinds(ri_loc) == istate_loc .and. self%jinds(ri_loc) == iproj_loc )then
                    fallback_ref_loc = ri_loc
                    exit
                endif
            enddo
            if( fallback_ref_loc == 0 )then
                do ri_loc = 1, self%nrefs
                    if( self%sinds(ri_loc) == istate_loc )then
                        fallback_ref_loc = ri_loc
                        exit
                    endif
                enddo
            endif
            if( fallback_ref_loc == 0 ) fallback_ref_loc = 1
            if( self%loc_tab(fallback_ref_loc,iptcl_loc)%inpl <= 0 )then
                fallback_state = self%sinds(fallback_ref_loc)
                fallback_proj  = self%jinds(fallback_ref_loc)
                fallback_ref_full = (fallback_state-1)*self%p_ptr%nspace + fallback_proj
                sh_seed = 0.
                if( self%p_ptr%l_doshift ) sh_seed = o_prev_loc%get_2Dshift()
                call self%b_ptr%pftc%gen_objfun_vals(fallback_ref_full, self%pinds(iptcl_loc), sh_seed, dists_inpl)
                dists_inpl = eulprob_dist_switch(dists_inpl, self%p_ptr%cc_objfun)
                irot_loc = self%b_ptr%pftc%get_roind(360.-o_prev_loc%e3get())
                if( self%p_ptr%l_prob_inpl )then
                    inpl_athres_state = calc_athres(self%b_ptr%spproj_field, 'dist_inpl', self%p_ptr%prob_athres, state=fallback_state)
                    irot_loc = angle_sampling(dists_inpl, dists_inpl_sorted, inds_sorted, inpl_athres_state, self%p_ptr%prob_athres)
                else
                    irot_loc = minloc(dists_inpl, dim=1)
                endif
                if( irot_loc < 1 .or. irot_loc > self%b_ptr%pftc%get_nrots() ) irot_loc = 1
                self%loc_tab(fallback_ref_loc,iptcl_loc)%dist   = dists_inpl(irot_loc)
                self%loc_tab(fallback_ref_loc,iptcl_loc)%inpl   = irot_loc
                self%loc_tab(fallback_ref_loc,iptcl_loc)%x      = sh_seed(1)
                self%loc_tab(fallback_ref_loc,iptcl_loc)%y      = sh_seed(2)
                self%loc_tab(fallback_ref_loc,iptcl_loc)%has_sh = self%p_ptr%l_doshift
            endif
            call mark_ref_touched(iptcl_loc, fallback_ref_loc)
            call o_prev_loc%kill
        end subroutine seed_fallback_if_empty

        subroutine advance_ref_head(iref)
            integer, intent(in) :: iref
            integer :: mloc, sloc, cand
            mloc = graph%ref_counts(iref)
            if( graph%ref_pos(iref) > mloc )then
                frontier%iref_dist(iref) = huge_val
                return
            endif
            sloc = graph%ref_offsets(iref)
            do while( graph%ref_pos(iref) <= mloc )
                cand = graph%ref_list(sloc + graph%ref_pos(iref) - 1)
                if( frontier%ptcl_avail(cand) )then
                    frontier%iref_dist(iref) = graph%ref_dists(sloc + graph%ref_pos(iref) - 1)
                    return
                endif
                graph%ref_pos(iref) = graph%ref_pos(iref) + 1
            enddo
            frontier%iref_dist(iref) = huge_val
        end subroutine advance_ref_head

        subroutine sync_frontier_ref(iref)
            integer, intent(in) :: iref
            if( frontier%sel_pos(iref) > 0 )then
                if( frontier%iref_dist(iref) >= huge_val )then
                    pos = frontier%sel_pos(iref)
                    if( pos < nsel )then
                        last_ref                   = frontier%sel_refs(nsel)
                        frontier%sel_refs(pos)     = last_ref
                        frontier%sel_dists(pos)    = frontier%sel_dists(nsel)
                        frontier%sel_pos(last_ref) = pos
                    endif
                    frontier%sel_pos(iref) = 0
                    nsel = nsel - 1
                else
                    frontier%sel_dists(frontier%sel_pos(iref)) = frontier%iref_dist(iref)
                endif
            else
                if( frontier%iref_dist(iref) < huge_val )then
                    nsel                     = nsel + 1
                    frontier%sel_refs(nsel)  = iref
                    frontier%sel_dists(nsel) = frontier%iref_dist(iref)
                    frontier%sel_pos(iref)   = nsel
                endif
            endif
        end subroutine sync_frontier_ref

    end subroutine ref_assign_neigh

    subroutine kill_neigh( self )
        class(eul_prob_tab_neigh), intent(inout) :: self
        if( allocated(self%eval_touched_refs)   ) deallocate(self%eval_touched_refs)
        if( allocated(self%eval_touched_counts) ) deallocate(self%eval_touched_counts)
        self%eval_max_touched = 0
        call self%eul_prob_tab%kill
    end subroutine kill_neigh

end module simple_eul_prob_tab_neigh
