program simple_test_uniform_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: FROM_IP=2, TO_IP=8, FROM_IR=3, TO_IR=5, NP=TO_IP-FROM_IP+1, NR=TO_IR-FROM_IR+1
real    :: ori_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), norm_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), tmp_sum
integer :: iptcl, iref, best_ip(FROM_IP:TO_IP), best_ir(FROM_IR:TO_IR), num, ind, from_ind, to_ind
num = int(NP/NR)
call random_number(ori_tab)
print *, 'Original table:'
do iptcl = FROM_IP,TO_IP
    print *, ori_tab(iptcl,:)
enddo
norm_tab = ori_tab
do iptcl = FROM_IP,TO_IP
    tmp_sum           = sum(norm_tab(iptcl,:))
    norm_tab(iptcl,:) = norm_tab(iptcl,:) / tmp_sum
enddo
print *, 'Normalized table:'
do iptcl = FROM_IP,TO_IP
    print *, norm_tab(iptcl,:)
enddo
best_ir  = (/(iref,  iref =1,NR)/)
best_ip  = (/(iptcl, iptcl=1,NP)/)
call cluster_sort(NP, NR, norm_tab, NP, NR, best_ip, best_ir)
best_ip  = best_ip + FROM_IP - 1
best_ir  = best_ir + FROM_IR - 1
print *, best_ir
print *, best_ip
print *, 'Sorted table:'
do iptcl = FROM_IP,TO_IP
    print *, norm_tab(best_ip(iptcl),best_ir)
enddo
print *, 'total prob of best ', num, 'entries:'
do iref = FROM_IR,TO_IR
    ind      = TO_IR-iref
    to_ind   = TO_IP-ind*num
    from_ind = to_ind-num+1
    print *, from_ind, to_ind
    print *, sum(norm_tab(best_ip(from_ind:to_ind), best_ir(iref)))
enddo
end program simple_test_uniform_cluster
