!@descr: project micrograph support, fixed-cardinality subset optimization and non-pruning selection masks
module simple_mic_subset
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_defs,         only: dp
use simple_sp_project,   only: sp_project
use simple_oris,         only: oris
use simple_ori,          only: ori
use simple_sym,          only: sym, is_valid_pointgroup
use simple_eo_partition, only: eo_partition_status, EO_COMPLETE
use simple_srch_sort_loc, only: hpsort
implicit none
private
public :: build_mic_support, select_mic_support, apply_mic_subset

contains

    subroutine build_mic_support(project, mode, score_key, state, nspace, pgrp, support, populations, &
                                 eligible, particle_mic, status, message, bin_eulers)
        class(sp_project), intent(inout) :: project
        character(len=*), intent(in) :: mode, score_key, pgrp
        integer, intent(in) :: state, nspace
        real(dp), allocatable, intent(out) :: support(:,:)
        integer, allocatable, intent(out) :: populations(:), particle_mic(:)
        logical, allocatable, intent(out) :: eligible(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        real, allocatable, optional, intent(out) :: bin_eulers(:,:)
        integer, allocatable :: mic2stk(:), stk2mic(:)
        type(oris) :: directions
        type(sym) :: symmetry
        type(ori) :: orientation
        character(len=512) :: eo_message
        integer :: nmics, nstks, np, i, m, s, b, half, ncells, physical_index, physical_count
        real :: euler(3), weight
        logical :: coverage
        status = 1
        message = ''
        coverage = mode == 'quality_coverage'
        if( .not.coverage .and. mode /= 'particles' )then
            message = 'subset_mode must be quality_coverage or particles'; return
        endif
        if( state < 1 .or. nspace < 1 )then
            message = 'state and nspace must be positive'; return
        endif
        if( coverage .and. mod(nspace,2) /= 0 )then
            message = 'quality_coverage requires an even nspace'; return
        endif
        if( score_key /= 'corr' .and. score_key /= 'corr_cart' .and. score_key /= 'none' )then
            message = 'subset_score must be corr, corr_cart or none'; return
        endif
        nmics = project%os_mic%get_noris()
        nstks = project%os_stk%get_noris()
        np = project%os_ptcl3D%get_noris()
        if( nmics < 1 .or. nstks < 1 .or. nmics < nstks .or. np < 1 )then
            message = 'project needs micrograph provenance, stacks and particles'; return
        endif
        if( project%os_ptcl2D%get_noris() /= np )then
            message = 'ptcl2D and ptcl3D row counts differ'; return
        endif
        call project%get_mic2stk_inds(mic2stk, stk2mic)
        do s = 1,nstks
            m = stk2mic(s)
            if( m < 1 .or. m > nmics )then
                message = 'stack has no valid micrograph mapping'; return
            endif
            if( project%os_mic%get_int(m,'nptcls') /= project%os_stk%get_int(s,'nptcls') )then
                message = 'micrograph and stack particle counts disagree'; return
            endif
        enddo
        ncells = 1
        if( coverage )then
            if( .not.is_valid_pointgroup(pgrp) )then
                message = 'invalid point-group symmetry'; return
            endif
            if( eo_partition_status(project, eo_message) /= EO_COMPLETE )then
                message = 'quality_coverage requires complete existing project half provenance'; return
            endif
            call symmetry%new(pgrp)
            call directions%new(nspace, is_ptcl=.false.)
            call symmetry%build_refspiral(directions)
            ncells = 2*directions%get_noris()
        endif
        if( present(bin_eulers) )then
            allocate(bin_eulers(3,max(1,ncells/2)),source=0.)
            if( coverage )then
                do b = 1,directions%get_noris()
                    bin_eulers(:,b) = directions%get_euler(b)
                enddo
            endif
        endif
        allocate(support(ncells,nmics), source=0._dp)
        allocate(populations(nmics), source=0)
        allocate(particle_mic(np), source=0)
        allocate(eligible(np), source=.false.)
        do i = 1,np
            if( project%os_ptcl2D%get_state(i) <= 0 ) cycle
            if( project%os_ptcl3D%get_state(i) /= state ) cycle
            s = project%os_ptcl3D%get_int(i,'stkind')
            if( s < 1 .or. s > nstks )then
                message = 'active particle has invalid stkind'; return
            endif
            if( project%os_ptcl2D%get_int(i,'stkind') /= s )then
                message = 'ptcl2D and ptcl3D stack mappings disagree'; return
            endif
            if( i < project%os_stk%get_fromp(s) .or. i > project%os_stk%get_top(s) )then
                message = 'particle lies outside its logical stack range'; return
            endif
            m = stk2mic(s)
            if( project%os_mic%get_state(m) <= 0 .or. project%os_stk%get_state(s) <= 0 ) cycle
            if( .not.project%os_stk%isthere(s,'nptcls_stk') )then
                message = 'eligible particle stack lacks its physical image count'; return
            endif
            physical_count = project%os_stk%get_int(s,'nptcls_stk')
            physical_index = project%os_ptcl3D%get_int(i,'indstk')
            if( physical_index < 1 .or. physical_index > physical_count )then
                message = 'eligible particle has no valid physical stack index'; return
            endif
            if( project%os_ptcl2D%get_int(i,'indstk') /= physical_index )then
                message = 'ptcl2D and ptcl3D physical image indices disagree'; return
            endif
            b = 1
            weight = 1.
            if( coverage )then
                if( project%os_ptcl3D%get_int(i,'updatecnt') <= 0 ) cycle
                euler = project%os_ptcl3D%get_euler(i)
                if( .not.all(ieee_is_finite(euler)) ) cycle
                if( euler(2) < 0. .or. euler(2) > 180. ) cycle
                if( score_key /= 'none' ) weight = project%os_ptcl3D%get(i,score_key)
                if( .not.ieee_is_finite(weight) ) cycle
                if( weight <= 0. ) cycle
                half = project%os_ptcl3D%get_eo(i)
                if( half < 0 .or. half > 1 )then
                    message = 'eligible particle has an invalid half label'; return
                endif
                if( project%os_ptcl2D%get_eo(i) /= half )then
                    message = 'eligible particle halves disagree between ptcl2D and ptcl3D'; return
                endif
                call project%os_ptcl3D%get_ori(i,orientation)
                b = 2*(symmetry%find_closest_proj(directions,orientation)-1) + half + 1
            endif
            eligible(i) = .true.
            particle_mic(i) = m
            populations(m) = populations(m) + 1
            support(b,m) = support(b,m) + real(weight,dp)
        enddo
        call orientation%kill
        call directions%kill
        call symmetry%kill
        status = 0
    end subroutine build_mic_support

    pure real(dp) function mic_subset_objective(total, scale) result(value)
        real(dp), intent(in) :: total(:), scale
        value = sum(log(1._dp + total/scale))
    end function mic_subset_objective

    subroutine select_mic_support(support, nselect, count_only, max_nodes, chosen, value, upper_bound, &
                                  proven, nodes, status)
        real(dp), intent(in) :: support(:,:)
        integer, intent(in) :: nselect, max_nodes
        logical, intent(in) :: count_only
        logical, allocatable, intent(out) :: chosen(:)
        real(dp), intent(out) :: value, upper_bound
        logical, intent(out) :: proven
        integer, intent(out) :: nodes, status
        integer, allocatable :: indices(:), selected(:), best(:)
        real(dp), allocatable :: h(:,:), suffix(:,:), largest(:,:), total(:), masses(:)
        real(dp) :: scale, gain, best_gain, tolerance
        integer :: n, nc, m, k, pick
        logical :: limited
        allocate(chosen(size(support,2)), source=.false.)
        value = 0._dp; upper_bound = 0._dp; proven = .false.; nodes = 0; status = 1
        if( nselect < 1 .or. max_nodes < 1 .or. size(support,1) < 1 ) return
        if( .not.all(ieee_is_finite(support)) ) return
        if( any(support < 0._dp) ) return
        n = count(sum(support,dim=1) > 0._dp)
        if( nselect > n ) return
        indices = pack([(m,m=1,size(support,2))], sum(support,dim=1) > 0._dp)
        masses = sum(support,dim=1)
        call hpsort(indices, preferred)
        if( count_only )then
            chosen(indices(1:nselect)) = .true.
            value = sum(masses(indices(1:nselect))); upper_bound = value; proven = .true.; status = 0
            return
        endif
        nc = size(support,1)
        h = support(:,indices)
        scale = max(tiny(1._dp), sum(h)*real(nselect,dp)/real(n,dp)/real(nc,dp))
        allocate(total(nc), source=0._dp)
        allocate(selected(nselect), best(nselect), source=0)
        do k = 1,nselect
            pick = 0; best_gain = -huge(1._dp)
            do m = 1,n
                if( any(best == m) ) cycle
                gain = mic_subset_objective(total+h(:,m),scale)
                if( gain > best_gain )then
                    best_gain = gain; pick = m
                endif
            enddo
            best(k) = pick; total = total+h(:,pick)
        enddo
        value = mic_subset_objective(total,scale)
        allocate(suffix(nc,n+1), largest(nc,n+1), source=0._dp)
        do m = n,1,-1
            suffix(:,m) = suffix(:,m+1)+h(:,m)
            largest(:,m) = max(largest(:,m+1),h(:,m))
        enddo
        upper_bound = mic_subset_objective(min(suffix(:,1),real(nselect,dp)*largest(:,1)),scale)
        tolerance = 1.e-12_dp*max(1._dp,upper_bound)
        limited = .false.
        total = 0._dp
        call search(1,0)
        chosen(indices(best)) = .true.
        proven = .not.limited
        if( proven ) upper_bound = value
        status = 0
      contains
        logical function preferred(a,b)
            integer, intent(in) :: a, b
            if( masses(a) == masses(b) )then
                preferred = a < b
            else
                preferred = masses(a) > masses(b)
            endif
        end function preferred

        recursive subroutine search(first, depth)
            integer, intent(in) :: first, depth
            real(dp) :: bound, score
            integer :: remaining, i
            remaining = nselect-depth
            if( n-first+1 < remaining ) return
            if( nodes >= max_nodes )then
                limited = .true.; return
            endif
            nodes = nodes+1
            bound = mic_subset_objective(total+min(suffix(:,first),real(remaining,dp)*largest(:,first)),scale)
            if( bound <= value+tolerance ) return
            if( remaining == 0 )then
                score = mic_subset_objective(total,scale)
                if( score > value+tolerance )then
                    value = score; best = selected
                endif
                return
            endif
            do i = first,n-remaining+1
                if( limited ) exit
                selected(depth+1) = i
                total = total+h(:,i)
                call search(i+1,depth+1)
                total = total-h(:,i)
            enddo
        end subroutine search
    end subroutine select_mic_support

    subroutine apply_mic_subset(project, chosen, eligible, particle_mic)
        class(sp_project), intent(inout) :: project
        logical, intent(in) :: chosen(:), eligible(:)
        integer, intent(in) :: particle_mic(:)
        integer, allocatable :: mic2stk(:), stk2mic(:)
        integer :: i, m
        logical :: keep
        call project%get_mic2stk_inds(mic2stk,stk2mic)
        do m = 1,size(chosen)
            if( .not.chosen(m) ) call project%os_mic%set_state(m,0)
        enddo
        do i = 1,size(stk2mic)
            if( .not.chosen(stk2mic(i)) ) call project%os_stk%set_state(i,0)
        enddo
        do i = 1,size(eligible)
            keep = .false.
            if( eligible(i) )then
                m = particle_mic(i)
                if( m >= 1 .and. m <= size(chosen) ) keep = chosen(m)
            endif
            if( keep ) cycle
            call project%os_ptcl2D%set_state(i,0)
            call project%os_ptcl3D%set_state(i,0)
        enddo
        call project%os_out%kill
        if( project%projinfo%isthere('sigma2_state') ) call project%projinfo%delete_entry(1,'sigma2_state')
    end subroutine apply_mic_subset

end module simple_mic_subset
