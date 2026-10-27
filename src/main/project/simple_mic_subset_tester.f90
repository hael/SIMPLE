!@descr: micrograph subset optimization against closed forms and an exhaustive oracle, with project mask invariants
module simple_mic_subset_tester
use simple_defs,          only: dp
use simple_test_utils,    only: assert_true, assert_int, assert_double
use simple_mic_subset,    only: select_mic_support, build_mic_support, rank_mic_support, apply_mic_subset
use simple_sp_project,    only: sp_project
use simple_eo_partition,  only: eo_transaction
implicit none
private
public :: run_all_mic_subset_tests

contains

    subroutine run_all_mic_subset_tests()
        call test_counts_and_refusals
        call test_global_search
        call test_exhaustive_oracle
        call test_project_mask
        call test_shared_coverage
    end subroutine run_all_mic_subset_tests

    subroutine test_shared_coverage()
        use simple_oris,        only: oris
        use simple_sym,         only: sym
        use simple_ptcl_subset, only: select_ptcl_subset
        class(sp_project), allocatable :: project
        class(oris), allocatable :: directions
        type(sym) :: symmetry
        real(dp), allocatable :: support(:,:), full_support(:,:)
        real(dp) :: expected(16,3)
        integer, allocatable :: populations(:), particle_mic(:)
        logical, allocatable :: eligible(:), selected(:)
        logical :: changed
        integer :: i, bin, cell, half, status, stack, target, best(16), halves(6)
        character(len=512) :: message
        allocate(project,directions)
        call fixture(project)
        call project%os_ptcl2D%set_all2single('state',1)
        call project%os_ptcl3D%set_all2single('state',2)
        call project%os_mic%set_all2single('state',1)
        call project%os_stk%set_all2single('state',1)
        call eo_transaction(project,changed)
        call directions%new(8,is_ptcl=.false.)
        call symmetry%new('c1')
        call symmetry%build_refspiral(directions)
        expected = 0._dp
        best = 0
        do i = 1,6
            bin = mod(i-1,2)+1
            call project%os_ptcl3D%set_euler(i,directions%get_euler(bin))
            call project%os_ptcl3D%set(i,'corr',real(i))
            half = project%os_ptcl3D%get_eo(i)
            halves(i) = half
            cell = 2*(bin-1)+half+1
            stack = project%os_ptcl3D%get_int(i,'stkind')
            expected(cell,stack) = expected(cell,stack) + real(i,dp)
            best(cell) = i
        enddo
        call build_mic_support(project,'quality_coverage','corr',2,8,'c1',support,populations, &
            eligible,particle_mic,status,message)
        call assert_int(0,status,'shared coverage fixture builds micrograph support')
        if( status == 0 )then
            call assert_true(all(abs(support-expected)<1.e-12_dp), &
                'micrograph support matches independently assigned directions, halves and scores')
        endif
        call select_ptcl_subset(project,4,8,'c1','corr',2,selected,status,message,rebalance_eo=.false.)
        call assert_int(0,status,'same project supports particle coverage selection')
        call assert_true(all(selected .eqv. [.false.,.false.,.true.,.true.,.true.,.true.]), &
            'shared classification preserves original top-two-per-direction particle selection')
        full_support = support
        call rank_mic_support(project,2,8,'c1','corr',2,populations,eligible,particle_mic, &
            support,target,status,message)
        call assert_int(0,status,'micrograph ranking calls the original particle engine in both halves')
        call assert_int(4,target,'six eligible particles times two of three micrographs gives target four')
        expected = 0._dp
        do cell = 1,16
            i = best(cell)
            if( i == 0 ) cycle
            stack = project%os_ptcl3D%get_int(i,'stkind')
            expected(cell,stack) = real(i,dp)
        enddo
        call assert_true(all(abs(support-expected)<1.e-12_dp), &
            'half-specific engine recommendations match the independent best-per-cell oracle')
        call assert_true(all(populations==2) .and. all(eligible),'recommendations preserve the full output eligibility pool')
        do i = 1,6
            call assert_int(halves(i),project%os_ptcl3D%get_eo(i),'engine calls preserve source 3D halves')
            call assert_int(halves(i),project%os_ptcl2D%get_eo(i),'engine calls preserve source 2D halves')
            call assert_int(2,project%os_ptcl3D%get_state(i),'engine calls preserve source state labels')
            call assert_double(real(i,dp),real(project%os_ptcl3D%get(i,'corr'),dp),'engine calls preserve source scores')
        enddo
        call rank_mic_support(project,3,8,'c1','corr',2,populations,eligible,particle_mic, &
            support,target,status,message)
        call assert_int(0,status,'selecting all candidate micrographs succeeds')
        call assert_int(6,target,'all micrographs recommend all eligible particles')
        call assert_true(all(abs(support-full_support)<1.e-12_dp),'full target retains original raw support')
        call rank_mic_support(project,1,8,'c1','none',2,populations,eligible,particle_mic, &
            support,target,status,message)
        call assert_int(0,status,'unweighted ranking adapts unit scores without modifying the particle engine')
        call assert_int(2,target,'small budget reserves one recommendation per populated half')
        call assert_double(1._dp,sum(support(1::2,:)),'one even-half unit recommendation')
        call assert_double(1._dp,sum(support(2::2,:)),'one odd-half unit recommendation')
        do i = 1,6
            call assert_double(real(i,dp),real(project%os_ptcl3D%get(i,'corr'),dp),'unit ranking leaves source scores intact')
        enddo
        do i = 1,6
            if( halves(i) == 0 ) cycle
            call project%os_ptcl2D%set_state(i,0)
            call project%os_ptcl3D%set_state(i,0)
        enddo
        call build_mic_support(project,'quality_coverage','none',2,8,'c1',support,populations, &
            eligible,particle_mic,status,message)
        call assert_int(0,status,'one populated half remains a valid candidate pool')
        call rank_mic_support(project,1,8,'c1','none',2,populations,eligible,particle_mic, &
            support,target,status,message)
        call assert_int(0,status,'half quota transfers to the available half')
        call assert_int(2,target,'one-half fixture derives two recommendations for one micrograph')
        call assert_double(2._dp,sum(support(1::2,:)),'available half receives the full target')
        call assert_double(0._dp,sum(support(2::2,:)),'engine never invents recommendations in an absent half')
        call rank_mic_support(project,4,8,'c1','corr',2,populations,eligible,particle_mic, &
            support,target,status,message)
        call assert_true(status/=0,'micrograph ranking refuses a budget larger than the eligible pool')
        call directions%kill
        call symmetry%kill
        call project%kill
    end subroutine test_shared_coverage

    subroutine test_counts_and_refusals()
        real(dp) :: h(1,4), value, upper
        logical, allocatable :: chosen(:)
        logical :: proven
        integer :: status, nodes
        h(1,:) = [3._dp,6._dp,6._dp,0._dp]
        call select_mic_support(h,2,.true.,1,chosen,value,upper,proven,nodes,status)
        call assert_int(0,status,'count-mode selection succeeds')
        call assert_true(all(chosen .eqv. [.false.,.true.,.true.,.false.]),'largest counts selected')
        call assert_double(12._dp,value,'closed-form maximum particle yield')
        call assert_true(proven .and. upper==value,'count optimum is certified')
        call select_mic_support(h,1,.true.,1,chosen,value,upper,proven,nodes,status)
        call assert_true(chosen(2) .and. .not.chosen(3),'count ties use original micrograph order')
        call select_mic_support(h,4,.true.,1,chosen,value,upper,proven,nodes,status)
        call assert_true(status/=0 .and. .not.any(chosen),'reject budget exceeding nonempty candidates')
        call select_mic_support(h,0,.false.,1,chosen,value,upper,proven,nodes,status)
        call assert_true(status/=0,'reject zero micrograph budget')
        h(1,1) = -1._dp
        call select_mic_support(h,1,.false.,1,chosen,value,upper,proven,nodes,status)
        call assert_true(status/=0,'reject negative support')
        h = 0._dp
        h(1,3) = 4._dp
        call select_mic_support(h,2,.false.,100,chosen,value,upper,proven,nodes,status, &
            candidate_mask=[.true.,.true.,.true.,.false.])
        call assert_int(0,status,'zero-recommendation micrographs can fill the exact budget')
        call assert_true(all(chosen .eqv. [.true.,.false.,.true.,.false.]), &
            'fillers use original micrograph order and never revive ineligible micrographs')
        call assert_true(proven,'zero-support filler objective is certified')
    end subroutine test_counts_and_refusals

    subroutine test_global_search()
        real(dp) :: h(6,3), value, upper, scale, expected
        logical, allocatable :: chosen(:)
        logical :: proven
        integer :: status, nodes
        h = 0._dp
        h(1:4,1) = .9_dp
        h([1,2,5],2) = 1._dp
        h([3,4,6],3) = 1._dp
        scale = (3.6_dp+3._dp+3._dp)*2._dp/3._dp/6._dp
        expected = 6._dp*log(1._dp+1._dp/scale)
        call select_mic_support(h,2,.false.,1000,chosen,value,upper,proven,nodes,status)
        call assert_int(0,status,'global search succeeds')
        call assert_true(all(chosen .eqv. [.false.,.true.,.true.]),'global search escapes greedy trap')
        call assert_double(expected,value,'complementary six-cell closed-form objective',1000000._dp)
        call assert_true(proven,'complete search certifies the optimum')
        call select_mic_support(h,2,.false.,1,chosen,value,upper,proven,nodes,status)
        call assert_true(status==0 .and. .not.proven,'limited search reports an unproven incumbent')
        call assert_int(2,count(chosen),'limited search still meets exact cardinality')
        call assert_true(upper>=expected .and. value<=expected+1.e-10_dp,'limited search retains a valid bound')
    end subroutine test_global_search

    subroutine test_exhaustive_oracle()
        real(dp) :: h(4,6), value, upper, scale, expected, score, total(4)
        logical, allocatable :: chosen(:)
        logical :: proven
        integer :: k, mask, m, status, nodes
        h = reshape([2._dp,0._dp,1._dp,0._dp, 0._dp,3._dp,0._dp,1._dp, &
                     1._dp,1._dp,1._dp,1._dp, 0._dp,0._dp,4._dp,2._dp, &
                     3._dp,0._dp,0._dp,2._dp, 0._dp,2._dp,2._dp,0._dp],shape(h))
        do k = 1,6
            scale = sum(h)*real(k,dp)/6._dp/4._dp
            expected = -1._dp
            do mask = 0,63
                if( popcnt(mask) /= k ) cycle
                total = 0._dp
                do m = 1,6
                    if( btest(mask,m-1) ) total = total+h(:,m)
                enddo
                score = sum(log(1._dp+total/scale))
                expected = max(expected,score)
            enddo
            call select_mic_support(h,k,.false.,10000,chosen,value,upper,proven,nodes,status)
            call assert_true(status==0 .and. proven,'exhaustive fixture finishes')
            call assert_int(k,count(chosen),'oracle fixture cardinality')
            call assert_double(expected,value,'search agrees with independent bitmask enumeration',1000000._dp)
        enddo
    end subroutine test_exhaustive_oracle

    subroutine fixture(project)
        class(sp_project), intent(inout) :: project
        integer :: i, s
        call project%os_mic%new(3,is_ptcl=.false.)
        call project%os_stk%new(3,is_ptcl=.false.)
        call project%os_ptcl2D%new(6,is_ptcl=.true.)
        call project%os_ptcl3D%new(6,is_ptcl=.true.)
        call project%os_out%new(1,is_ptcl=.false.)
        call project%os_out%set(1,'vol','/parent/volume.mrc')
        call project%projinfo%new(1,is_ptcl=.false.)
        call project%projinfo%set(1,'sigma2_state','/parent/sigma.state')
        do s = 1,3
            call project%os_mic%set(s,'state',1)
            call project%os_mic%set(s,'nptcls',2)
            call project%os_stk%set(s,'state',1)
            call project%os_stk%set(s,'nptcls',2)
            call project%os_stk%set(s,'nptcls_stk',2)
            call project%os_stk%set(s,'fromp',2*s-1)
            call project%os_stk%set(s,'top',2*s)
        enddo
        do i = 1,6
            s = (i+1)/2
            call project%os_ptcl2D%set_state(i,1)
            call project%os_ptcl3D%set_state(i,2)
            call project%os_ptcl2D%set(i,'stkind',s)
            call project%os_ptcl3D%set(i,'stkind',s)
            call project%os_ptcl2D%set(i,'indstk',mod(i-1,2)+1)
            call project%os_ptcl3D%set(i,'indstk',mod(i-1,2)+1)
            call project%os_ptcl3D%set(i,'pind',100+i)
            call project%os_ptcl3D%set(i,'updatecnt',1)
            call project%os_ptcl3D%set(i,'corr',.5)
            call project%os_ptcl3D%set_euler(i,[real(15*i),30.,10.])
        enddo
        call project%os_ptcl2D%set_state(2,0)
        call project%os_ptcl3D%set_state(2,0)
        call project%os_ptcl3D%set_state(4,3)
        call project%os_mic%set_state(3,0)
        call project%os_stk%set_state(3,0)
    end subroutine fixture

    subroutine test_project_mask()
        class(sp_project), allocatable :: project
        real(dp), allocatable :: support(:,:)
        integer, allocatable :: populations(:), particle_mic(:)
        logical, allocatable :: eligible(:), chosen(:)
        real(dp) :: value, upper
        real :: euler(3)
        integer :: status, nodes, half
        logical :: proven, changed
        character(len=512) :: message
        allocate(project)
        call fixture(project)
        call project%os_ptcl3D%set_all2single('updatecnt',0)
        call build_mic_support(project,'particles','corr',2,8,'c1',support,populations, &
                               eligible,particle_mic,status,message)
        call assert_int(0,status,'count mode needs neither poses nor half provenance')
        call assert_true(all(populations==[1,1,0]),'inactive particles, other states and rejected mics excluded')
        call assert_int(3,project%os_ptcl3D%get_state(4),'building support leaves source states unchanged')
        call assert_int(1,project%os_out%get_noris(),'building support leaves inherited outputs unchanged')
        call build_mic_support(project,'quality_coverage','corr',2,8,'c1',support,populations, &
                               eligible,particle_mic,status,message)
        call assert_true(status/=0,'coverage refuses absent half provenance')
        call project%os_ptcl3D%set_all2single('updatecnt',1)
        call eo_transaction(project,changed)
        call project%os_ptcl2D%set(1,'stkind',2)
        call build_mic_support(project,'particles','none',2,8,'c1',support,populations, &
                               eligible,particle_mic,status,message)
        call assert_true(status/=0,'inconsistent particle stack provenance is refused')
        call project%os_ptcl2D%set(1,'stkind',1)
        call build_mic_support(project,'quality_coverage','corr',2,8,'c1',support,populations, &
                               eligible,particle_mic,status,message)
        call assert_int(0,status,'coverage consumes existing half provenance')
        if( status/=0 )then
            call project%kill
            return
        endif
        call assert_double(1._dp,sum(support),'two eligible particles contribute their stored score')
        call select_mic_support(support,1,.false.,1000,chosen,value,upper,proven,nodes,status)
        chosen = [.true.,.false.,.false.]
        euler = project%os_ptcl3D%get_euler(1)
        half = project%os_ptcl3D%get_eo(1)
        call apply_mic_subset(project,chosen,eligible,particle_mic)
        call assert_int(6,project%os_ptcl3D%get_noris(),'mask never renumbers particles')
        call assert_int(3,project%os_mic%get_noris(),'mask never compacts micrograph rows')
        call assert_int(2,project%os_ptcl3D%get_state(1),'retained structural label preserved')
        call assert_int(0,project%os_ptcl3D%get_state(2),'rejected particle stays rejected')
        call assert_int(0,project%os_ptcl3D%get_state(3),'unselected micrograph particles deselected')
        call assert_int(1,project%os_ptcl3D%get_int(1,'indstk'),'physical stack index preserved')
        call assert_int(101,project%os_ptcl3D%get_int(1,'pind'),'original particle identity preserved')
        call assert_int(half,project%os_ptcl3D%get_eo(1),'particle half preserved')
        call assert_true(all(project%os_ptcl3D%get_euler(1)==euler),'pose preserved')
        call assert_int(0,project%os_out%get_noris(),'subset clears parent reconstruction and weight outputs')
        call assert_true(.not.project%projinfo%isthere('sigma2_state'),'subset clears parent sigma registration')
        call project%os_mic%kill
        call build_mic_support(project,'particles','none',2,8,'c1',support,populations, &
                               eligible,particle_mic,status,message)
        call assert_true(status/=0,'missing micrograph provenance is refused')
        call project%kill
    end subroutine test_project_mask

end module simple_mic_subset_tester
