program simple_test_uniform_sort
include 'simple_lib.f08'
implicit none
integer, parameter :: FROM_IP=2, TO_IP=8, FROM_IR=1, TO_IR=3, NP=TO_IP-FROM_IP+1, NR=TO_IR-FROM_IR+1
real    :: ori_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), norm_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), tmp_sum
integer :: iptcl, iref, best_ip(FROM_IP:TO_IP), best_ir(FROM_IP:TO_IP), num, ind, from_ind, to_ind
logical :: mask_ir(FROM_IR:TO_IR)
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
mask_ir  = .false.
best_ip  = (/(iptcl, iptcl=1,NP)/)
best_ir  = FROM_IR
call uniform_sort(NP, NR, norm_tab, NP, best_ip, best_ir, mask_ir)
best_ip  = best_ip + FROM_IP - 1
best_ir  = best_ir + FROM_IR - 1
print *, best_ip
print *, best_ir
print *, 'Sorted table:'
do iptcl = FROM_IP,TO_IP
    print *, norm_tab(best_ip(iptcl),:)
enddo
end program simple_test_uniform_sort
