program simple_test_uniform_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: NP = 13, NR = 4
real    :: ori_tab(NP, NR), norm_tab(NP, NR), tmp_sum
integer :: iptcl, iref, best_id(NP), best_ir(NR)
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
call cluster_sort(norm_tab, NP, NR, best_id, best_ir)
print *, best_ir
print *, best_id
print *, 'Sorted table:'
do iptcl = 1,NP
    print *, norm_tab(best_id(iptcl),best_ir)
enddo
contains
    subroutine cluster_sort( tab, to_ii, to_ir, cur_id, cur_ir )
        real,    intent(in)    :: tab(NP, NR)
        integer, intent(in)    :: to_ii
        integer, intent(in)    :: to_ir
        integer, intent(inout) :: cur_id(NP)
        integer, intent(inout) :: cur_ir(NR)
        integer :: ir, ip, tmp_id(to_ii), tmp_i, orig_id(to_ii), ir_best
        real    :: best_sum, sum_prob
        if( to_ii <= int(NP/NR) ) return
        best_sum = 0.
        orig_id  = cur_id(1:to_ii)
        ir_best  = to_ir
        do ir = 1,to_ir
            tmp_id = (/(ip, ip=1,to_ii)/)
            call hpsort_ind(tmp_id, tab(orig_id,cur_ir(ir)))
            sum_prob = sum(tab(orig_id(tmp_id(to_ii-int(NP/NR)+1:to_ii)),cur_ir(ir)))
            if( sum_prob > best_sum )then
                best_sum        = sum_prob
                cur_id(1:to_ii) = orig_id(tmp_id)
                ir_best         = ir
            endif
        enddo
        ! swapping the last with the current best ir
        tmp_i           = cur_ir(to_ir)
        cur_ir(to_ir)   = cur_ir(ir_best)
        cur_ir(ir_best) = tmp_i
        call cluster_sort( tab, to_ii - int(NP/NR), to_ir - 1, cur_id, cur_ir )
    end subroutine cluster_sort

    ! sorting rarr, but only keep the sorted indeces
    subroutine hpsort_ind( iarr, rarr )
        integer, intent(inout) :: iarr(:)
        real,    intent(in)    :: rarr(:)
        integer :: i, ir, j, l, ia, n
        real    :: ra
        n = size(rarr)
        if( n /= size(iarr) )&
        &call simple_exception('nonconforming array sizes; hpsort_6', __FILENAME__ , __LINE__)
        if( n < 2 ) return
        l  = n/2+1
        ir = n
        do
            if(l > 1)then
                l  = l-1
                ia = iarr(l)
                ra = rarr(ia)
            else
                ia = iarr(ir)
                ra = rarr(ia)
                iarr(ir) = iarr(1)
                ir = ir-1
                if(ir == 1)then
                    iarr(1) = ia
                    return
                endif
            endif
            i = l
            j = l+l
            do while(j <= ir)
                if(j < ir) then
                    if(rarr(iarr(j)) < rarr(iarr(j+1))) j = j+1
                endif
                if(ra < rarr(iarr(j)))then
                    iarr(i) = iarr(j)
                    i = j
                    j = j+j
                else
                    j = ir+1
                endif
                iarr(i) = ia
            end do
        end do
    end subroutine hpsort_ind
end program simple_test_uniform_cluster
