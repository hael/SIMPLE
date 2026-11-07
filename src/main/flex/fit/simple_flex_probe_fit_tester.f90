!@descr: white-box test of the flex probe-part codec: one fit written by the workers' part writer folds back unchanged
submodule (simple_flex_probe_fit) simple_flex_probe_fit_tester
use simple_test_utils,      only: assert_true
use simple_syslib,          only: del_file
use simple_rnd,             only: seed_rnd_fixed
use simple_flex_pca_rounds, only: flex_pca_rounds_shmem
implicit none

contains

    module subroutine run_all_flex_probe_fit_tests()
        write(*,'(A)') '**** running all flex probe-fit tests ****'
        call test_probe_part_codec()
    end subroutine run_all_flex_probe_fit_tests

    !> every payload block of one fit survives the writer and the reduce bit for bit
    subroutine test_probe_part_codec()
        logical :: pass
        type(flex_probe_part) :: src, dst
        type(flex_pca_rounds_shmem) :: rounds
        class(parameters), allocatable :: params
        type(string) :: fname, tmp_fname
        integer :: i, j, k, q, p, funit
        real :: xr, xi
        src%ncomp  = 2
        src%nval   = 17
        src%nll_sum = 1.25d0
        allocate(src%cmat_e(3,2,2,2), src%cmat_o(3,2,2,2))
        allocate(src%rho_ex(3,2,2,2), src%rho_ox(3,2,2,2))
        allocate(src%rho_e(3,3,2,2), src%rho_o(3,3,2,2), src%gam_sum(2))
        allocate(src%mix_sr(2), src%mix_sm(2,2), src%mix_smm(2,2,2), src%mix_sainv(2,2))
        allocate(src%z_sub(5,2)); src%nz = 5
        allocate(src%kpk_e(3,4), src%kpk_o(3,4), src%rpk_e(2,4), src%rpk_o(2,4))
        dst%ncomp = src%ncomp
        allocate(dst%cmat_e(3,2,2,2), dst%cmat_o(3,2,2,2), source=cmplx(0.,0.))
        allocate(dst%rho_ex(3,2,2,2), dst%rho_ox(3,2,2,2), source=0.)
        allocate(dst%rho_e(3,3,2,2), dst%rho_o(3,3,2,2), source=0.)
        allocate(dst%gam_sum(2), source=0.d0)
        allocate(dst%mix_sr(2), dst%mix_sm(2,2), dst%mix_smm(2,2,2), dst%mix_sainv(2,2), source=0.d0)
        allocate(dst%z_sub(5,2), source=0.d0)
        allocate(dst%kpk_e(3,4), dst%kpk_o(3,4), source=0.)
        allocate(dst%rpk_e(2,4), dst%rpk_o(2,4), source=cmplx(0.,0.))
        call seed_rnd_fixed(20261008)
        do q = 1, 2
            do k = 1, 2
                do j = 1, 2
                    do i = 1, 3
                        call random_number(xr); call random_number(xi)
                        src%cmat_e(i,j,k,q) = cmplx(xr,xi)
                        call random_number(xr); call random_number(xi)
                        src%cmat_o(i,j,k,q) = cmplx(xr,xi)
                        call random_number(src%rho_ex(i,j,k,q))
                        call random_number(src%rho_ox(i,j,k,q))
                    end do
                end do
            end do
        end do
        call random_number(src%rho_e); call random_number(src%rho_o)
        call random_number(src%gam_sum); call random_number(src%mix_sr)
        call random_number(src%mix_sm); call random_number(src%mix_smm)
        call random_number(src%mix_sainv); call random_number(src%z_sub)
        call random_number(src%kpk_e); call random_number(src%kpk_o)
        do p = 1, 4
            do q = 1, 2
                call random_number(xr); call random_number(xi)
                src%rpk_e(q,p) = cmplx(xr,xi)
                call random_number(xr); call random_number(xi)
                src%rpk_o(q,p) = cmplx(xr,xi)
            end do
        end do
        allocate(params)
        params%numlen = 4
        fname = rounds%part_fname('probe', 1, params%numlen)
        call open_probe_part_write(fname, 1, funit, tmp_fname)
        call write_probe_part_fit(funit, src)
        call close_probe_part_write(funit, tmp_fname, fname)
        call reduce_probe_parts(params, rounds, dst)
        ! The reducer owns normal cleanup; keep the fixture idempotent if that lifecycle
        ! changes or an implementation returns after a successful read without deleting.
        call del_file(fname)
        pass = src%ncomp == dst%ncomp .and. src%nval == dst%nval .and. src%nll_sum == dst%nll_sum
        pass = pass .and. src%nz == dst%nz
        pass = pass .and. all(src%cmat_e == dst%cmat_e) .and. all(src%cmat_o == dst%cmat_o)
        pass = pass .and. all(src%rho_ex == dst%rho_ex) .and. all(src%rho_ox == dst%rho_ox)
        pass = pass .and. all(src%rho_e == dst%rho_e) .and. all(src%rho_o == dst%rho_o)
        pass = pass .and. all(src%gam_sum == dst%gam_sum)
        pass = pass .and. all(src%mix_sr == dst%mix_sr) .and. all(src%mix_sm == dst%mix_sm)
        pass = pass .and. all(src%mix_smm == dst%mix_smm) .and. all(src%mix_sainv == dst%mix_sainv)
        pass = pass .and. all(src%z_sub == dst%z_sub)
        pass = pass .and. all(src%kpk_e == dst%kpk_e) .and. all(src%kpk_o == dst%kpk_o)
        pass = pass .and. all(src%rpk_e == dst%rpk_e) .and. all(src%rpk_o == dst%rpk_o)
        call fname%kill
        call src%kill
        call dst%kill
        deallocate(params)
        call assert_true(pass, 'probe-part codec round trip is bit exact for every payload block')
    end subroutine test_probe_part_codec

end submodule simple_flex_probe_fit_tester
