!@descr: 3D strategy for coarse subspace search followed by corr-guided probabilistic block-tree descent
!        Phase 1 mirrors the coarse subspace search of strategy3D_greedy_sub to identify the peak trees.
!        Phase 2 descends each peak's tree probabilistically: at every internal node both children are
!        evaluated via pftc%gen_objfun_vals; the in-plane rotation is selected with angle_sampling
!        (prob-weighted, matching l_prob_inpl logic); the branching decision is drawn proportional to
!        exp(best_child_corr), i.e. the child with higher correlation is preferred stochastically.
module simple_strategy3D_ptree
use simple_core_module_api
use simple_strategy3D_alloc
use simple_strategy3D_utils
use simple_parameters,      only: parameters
use simple_oris,            only: oris
use simple_strategy3D,      only: strategy3D
use simple_strategy3D_srch, only: strategy3D_spec
implicit none

public :: strategy3D_ptree
private
#include "simple_local_flags.inc"

type, extends(strategy3D) :: strategy3D_ptree
contains
    procedure :: new         => new_ptree
    procedure :: srch        => srch_ptree
    procedure :: kill        => kill_ptree
    procedure :: oris_assign => oris_assign_ptree
end type strategy3D_ptree

contains

    subroutine new_ptree( self, params, spec, build )
        use simple_builder, only: builder
        class(strategy3D_ptree), intent(inout) :: self
        class(parameters),           intent(in)    :: params
        class(strategy3D_spec),      intent(inout) :: spec
        class(builder),              intent(in)    :: build
        call self%s%new(params, spec, build)
        self%spec = spec
    end subroutine new_ptree

    subroutine srch_ptree( self, os, ithr )
        use simple_eul_prob_tab, only: angle_sampling, eulprob_dist_switch
        use simple_binary_tree,  only: bt_node
        class(strategy3D_ptree), intent(inout) :: self
        class(oris),                 intent(inout) :: os
        integer,                     intent(in)    :: ithr
        ! coarse search locals
        integer :: iref, isample, loc(1), iproj, ipeak
        real    :: inpl_corrs(self%s%nrots)
        ! tree descent locals
        integer              :: ntrees, itree, inode, iref_L, iref_R, nrefs_tree, istate
        integer              :: peak_refs(self%s%npeaks)
        logical, allocatable :: trees_done(:)
        real                 :: corrs_L(self%s%nrots), corrs_R(self%s%nrots)
        real                 :: sorted_L(self%s%nrots), sorted_R(self%s%nrots)
        integer              :: inds_L(self%s%nrots),  inds_R(self%s%nrots)
        integer              :: loc_L(1), loc_R(1)
        real                 :: best_corr_L, best_corr_R, corr_tmp, p_left, p_right
        type(bt_node)        :: node_cur, node_L, node_R, node_root
        if( os%get_state(self%s%iptcl) > 0 )then
            ! set thread index
            self%s%ithr = ithr
            ! prep
            call self%s%prep4srch
            ! shift search on previous best reference
            call self%s%inpl_srch_first
            ! ------------------------------------------------------------------
            ! Phase 1: coarse search over the subspace (mirrors greedy_sub)
            !          uses maxloc for the in-plane selection, consistent with
            !          the greedy coarse pass used to seed the neighborhoods.
            ! ------------------------------------------------------------------
            do isample = 1, self%s%nrefs_sub
                iref = s3D%srch_order_sub(isample, self%s%ithr)
                if( s3D%state_exists(s3D%proj_space_state(iref)) )then
                    if( self%s%p_ptr%l_sh_first )then
                        call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, self%s%xy_first, inpl_corrs)
                    else
                        call self%s%b_ptr%pftc%gen_objfun_vals(iref, self%s%iptcl, [0.,0.],         inpl_corrs)
                    endif
                    loc = maxloc(inpl_corrs)
                    call self%s%store_solution(iref, loc(1), inpl_corrs(loc(1)))
                endif
            end do
            ! identify the top-npeaks refs from the coarse search to select trees
            peak_refs = maxnloc(s3D%proj_space_corrs(:, self%s%ithr), self%s%npeaks)
            ! ------------------------------------------------------------------
            ! Phase 2: corr-guided probabilistic tree descent
            !          one descent per unique tree that contains a coarse peak.
            !          At each internal node, both children are evaluated with
            !          gen_objfun_vals + angle_sampling; the branch is drawn
            !          proportional to exp(best_child_corr).
            ! ------------------------------------------------------------------
            ntrees = self%s%b_ptr%block_tree%get_n_trees()
            allocate(trees_done(ntrees), source=.false.)
            nrefs_tree = 0
            do ipeak = 1, self%s%npeaks
                iproj = s3D%proj_space_proj(peak_refs(ipeak))
                if( .not. allocated(self%s%b_ptr%subspace_full2sub_map) )&
                &THROW_HARD('Probabilistic tree search requires subspace_full2sub_map to identify trees. Check builder construction.')
                itree = self%s%b_ptr%subspace_full2sub_map(iproj)
                if( itree < 1 .or. itree > ntrees )&
                &THROW_HARD('Invalid tree index mapped from peak reference. Check builder construction.')
                if( trees_done(itree) ) cycle
                trees_done(itree) = .true.
                ! start descent from the root of this tree
                node_root = self%s%b_ptr%block_tree%get_root_node(itree)
                inode     = node_root%node_idx
                do
                    if( inode == 0 ) exit
                    if( self%s%b_ptr%block_tree%is_leaf(itree, inode) ) exit
                    node_cur = self%s%b_ptr%block_tree%get_node(itree, inode)
                    if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
                    ! --- evaluate left child across all states ---
                    best_corr_L = -huge(1.0)
                    if( node_cur%left_idx /= 0 )then
                        node_L = self%s%b_ptr%block_tree%get_node(itree, node_cur%left_idx)
                        do istate = 1, self%s%nstates
                            iref_L = (istate - 1) * self%s%p_ptr%nspace + node_L%ref_idx
                            if( .not. s3D%state_exists(s3D%proj_space_state(iref_L)) ) cycle
                            if( self%s%p_ptr%l_sh_first )then
                                call self%s%b_ptr%pftc%gen_objfun_vals(iref_L, self%s%iptcl, self%s%xy_first, corrs_L)
                            else
                                call self%s%b_ptr%pftc%gen_objfun_vals(iref_L, self%s%iptcl, [0.,0.],         corrs_L)
                            endif
                            loc_L    = angle_sampling(eulprob_dist_switch(corrs_L, self%s%p_ptr%cc_objfun),&
                                &sorted_L, inds_L,&
                                &s3D%smpl_inpl_athres(s3D%proj_space_state(iref_L)), self%s%p_ptr%prob_athres)
                            corr_tmp = corrs_L(loc_L(1))
                            call self%s%store_solution(iref_L, loc_L(1), corr_tmp)
                            nrefs_tree  = nrefs_tree + 1
                            best_corr_L = max(best_corr_L, corr_tmp)
                        end do
                    endif
                    ! --- evaluate right child across all states ---
                    best_corr_R = -huge(1.0)
                    if( node_cur%right_idx /= 0 )then
                        node_R = self%s%b_ptr%block_tree%get_node(itree, node_cur%right_idx)
                        do istate = 1, self%s%nstates
                            iref_R = (istate - 1) * self%s%p_ptr%nspace + node_R%ref_idx
                            if( .not. s3D%state_exists(s3D%proj_space_state(iref_R)) ) cycle
                            if( self%s%p_ptr%l_sh_first )then
                                call self%s%b_ptr%pftc%gen_objfun_vals(iref_R, self%s%iptcl, self%s%xy_first, corrs_R)
                            else
                                call self%s%b_ptr%pftc%gen_objfun_vals(iref_R, self%s%iptcl, [0.,0.],         corrs_R)
                            endif
                            loc_R    = angle_sampling(eulprob_dist_switch(corrs_R, self%s%p_ptr%cc_objfun),&
                                &sorted_R, inds_R,&
                                &s3D%smpl_inpl_athres(s3D%proj_space_state(iref_R)), self%s%p_ptr%prob_athres)
                            corr_tmp = corrs_R(loc_R(1))
                            call self%s%store_solution(iref_R, loc_R(1), corr_tmp)
                            nrefs_tree  = nrefs_tree + 1
                            best_corr_R = max(best_corr_R, corr_tmp)
                        end do
                    endif
                    ! --- probabilistic branching decision ---
                    !     weights: exp(best_corr); child absent => weight ≈ 0
                    p_left  = exp(best_corr_L)
                    p_right = exp(best_corr_R)
                    if( sample_two(p_left, p_right) == 1 .and. node_cur%left_idx /= 0 )then
                        inode = node_cur%left_idx
                    else if( node_cur%right_idx /= 0 )then
                        inode = node_cur%right_idx
                    else
                        exit
                    endif
                end do
            end do
            if( allocated(trees_done) ) deallocate(trees_done)
            ! count evaluations (coarse refs + tree nodes visited)
            self%s%nrefs_eval = self%s%nrefs_sub + nrefs_tree
            ! pick the global best npeaks from all solutions accumulated above
            call extract_peak_oris(self%s)
            ! refine in-plane rotation at each peak orientation
            call self%s%inpl_srch_peaks
            ! finalise the orientation assignment
            call self%oris_assign
        else
            call os%reject(self%s%iptcl)
        endif

        contains

            !> Draw from {1,2} with probabilities proportional to (p1, p2).
            !  Falls back to uniform if both weights are zero.
            function sample_two( p1, p2 ) result(which)
                real, intent(in) :: p1, p2
                integer          :: which
                real             :: r, psum
                psum = p1 + p2
                if( psum <= 0.0 )then
                    which = merge(1, 2, ran3() < 0.5)
                    return
                endif
                r     = ran3()    ! uniform [0,1)
                which = merge(1, 2, r < p1 / psum)
            end function sample_two

    end subroutine srch_ptree

    subroutine oris_assign_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_ptree

    subroutine kill_ptree( self )
        class(strategy3D_ptree), intent(inout) :: self
        call self%s%kill
    end subroutine kill_ptree

end module simple_strategy3D_ptree
