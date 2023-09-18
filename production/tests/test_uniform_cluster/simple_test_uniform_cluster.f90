program simple_test_uniform_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: NP = 6, NR = 3
real    :: ori_tab(NP, NR), sorted_tab(NP, NR), sum_prob, best_sum, cur_col(NP)
integer :: i, j, best_id(NP), best_ir, cur_id(NP)
call random_number(ori_tab)
print *, 'Original table:'
do i = 1,NP
    print *, ori_tab(i,:)
enddo
sorted_tab = ori_tab
do i = 1,NP
    sum_prob = sum(sorted_tab(i,:))
    sorted_tab(i,:) = sorted_tab(i,:) / sum_prob
enddo
print *, 'Normalized table:'
do i = 1,NP
    print *, sorted_tab(i,:)
enddo
best_sum = 0.
do j = 1,NR
    cur_id  = (/(i,i=1,NP)/)
    cur_col = sorted_tab(:,j)
    call hpsort(cur_col, cur_id)
    call reverse(cur_id)
    call reverse(cur_col)
    sum_prob = sum(cur_col(1:int(NP/NR)))
    if( sum_prob > best_sum )then
        best_sum = sum_prob
        best_id  = cur_id
        best_ir  = j
    endif
enddo
print *, 'Sum of each column:'
print *, sum(sorted_tab, dim=1)
print *, best_sum
print *, best_ir
print *, best_id
end program simple_test_uniform_cluster
