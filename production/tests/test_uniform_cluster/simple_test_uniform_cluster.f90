program simple_test_uniform_cluster
include 'simple_lib.f08'
implicit none
integer, parameter :: FROM_IP=2, TO_IP=8, FROM_IR=1, TO_IR=3, NP=TO_IP-FROM_IP+1, NR=TO_IR-FROM_IR+1
real    :: ori_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), norm_tab(FROM_IP:TO_IP, FROM_IR:TO_IR), tmp_sum
integer :: iptcl, best_ir(FROM_IP:TO_IP)
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
best_ir = 1
call reg_uniform_sort(NR, best_ir)
do iptcl = FROM_IP, TO_IP
    print *, iptcl, ' -> ', best_ir(iptcl)
enddo
contains
    subroutine reg_uniform_sort( ncols, cur_ir )
        integer, intent(in)    :: ncols
        integer, intent(inout) :: cur_ir(FROM_IP:TO_IP)
        logical :: mask_ir(ncols), mask_ip(FROM_IP:TO_IP)
        integer :: ir, ip, max_ind_ir, max_ind_ip, max_ip(ncols)
        real    :: max_ir(ncols)
        mask_ir = .false.
        mask_ip = .true.
        do
            if( .not.(any(mask_ip)) ) return
            if( .not.(any(mask_ir)) ) mask_ir = .true.
            max_ir = 0.
            !$omp parallel do default(shared) proc_bind(close) schedule(static) private(ir,ip)
            do ir = 1, ncols
                if( mask_ir(ir) )then
                    do ip = FROM_IP, TO_IP
                        if( mask_ip(ip) .and. norm_tab(ip, ir) > max_ir(ir) )then
                            max_ir(ir) = norm_tab(ip, ir)
                            max_ip(ir) = ip
                        endif
                    enddo
                endif
            enddo
            !$omp end parallel do
            max_ind_ir = maxloc(max_ir, dim=1, mask=mask_ir)
            max_ind_ip = max_ip(max_ind_ir)
            cur_ir( max_ind_ip) = max_ind_ir
            mask_ip(max_ind_ip) = .false.
            mask_ir(max_ind_ir) = .false.
            print *, max_ind_ip, ' ----> ', max_ind_ir
        enddo
    end subroutine reg_uniform_sort
end program simple_test_uniform_cluster
