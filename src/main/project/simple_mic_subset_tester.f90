!@descr: micrograph subset optimization against closed forms and an exhaustive oracle, with project mask invariants
module simple_mic_subset_tester
use simple_defs,          only: dp
use simple_test_utils,    only: assert_true, assert_int, assert_double
use simple_mic_subset,    only: select_mic_support, build_mic_support, apply_mic_subset
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
    end subroutine run_all_mic_subset_tests

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
