program simple_test_uniform_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: NP = 13, NR = 4
real    :: ori_tab(NP, NR), norm_tab(NP, NR), tmp_sum
integer :: iptcl, iref, best_id(NP), best_ir(NR), cnt
call random_number(ori_tab)
print *, 'Original table:'
do iptcl = 1,NP
    print *, ori_tab(iptcl,:)
enddo
norm_tab = ori_tab
do iptcl = 1,NP
    tmp_sum           = sum(norm_tab(iptcl,:))
    norm_tab(iptcl,:) = norm_tab(iptcl,:) / tmp_sum
enddo
print *, 'Normalized table:'
do iptcl = 1,NP
    print *, norm_tab(iptcl,:)
enddo
best_ir  = (/(iref,  iref =1,NR)/)
best_id  = (/(iptcl, iptcl=1,NP)/)
call cluster_sort(NP, NR, norm_tab, NP, NR, best_id, best_ir)
print *, best_ir
print *, best_id
print *, 'Sorted table:'
do iptcl = 1,NP
    print *, norm_tab(best_id(iptcl),best_ir)
enddo
print *, 'total prob of best ', int(NP/NR), 'entries:'
cnt = NR-1
do iref = 1,NR
    print *, sum(norm_tab(best_id(NP-(cnt+1)*int(NP/NR)+1:NP-cnt*int(NP/NR)), best_ir(iref)))
    cnt = cnt - 1
enddo
end program simple_test_uniform_cluster
