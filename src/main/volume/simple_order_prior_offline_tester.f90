!@descr: unit tests for the offline order-prior harness: surrogate operator against its closed form, symmetry, micelle retention, threads
submodule (simple_order_prior_offline) simple_order_prior_offline_tester
!$ use omp_lib,      only: omp_get_max_threads, omp_set_num_threads
use simple_test_utils, only: assert_true, set_fixed_seed
implicit none

integer, parameter :: TBOX     = 48
integer, parameter :: TSEED    = 20261107
real,    parameter :: TSMPD    = 2.0
real,    parameter :: TMSKDIAM = 80.      !< a 20-voxel sphere radius
real,    parameter :: TLAMBDA  = 2.0
real,    parameter :: SYM_TOL  = 1.e-5    !< |<Au,v> - <u,Av>| / (|Au||v|): single-precision transform round-off
real,    parameter :: TEAM_TOL = 1.e-6    !< elementwise loops and fixed FFTW plans: equal to round-off
real,    parameter :: PROT_R   = 6.5      !< protein radius of the fixture (voxels)
real,    parameter :: MIC_R    = 12.      !< micelle shell radius and width (voxels)
real,    parameter :: MIC_W    = 2.

contains

    module subroutine run_all_order_prior_offline_tests()
        write(logfhandle,'(A)') '**** running all order prior offline tests ****'
        call test_constant_precision_is_closed_form()
        call test_operator_symmetry()
        call test_micelle_retention()
        call test_team_matches_serial()
    end subroutine run_all_order_prior_offline_tests

    !> with D constant the operator is pointwise, so its solution is the closed form h/(h + lambda_s (1-w));
    !! a CG stop at relative residual r bounds the relative error by cond * r, cond = 1 + lambda_rel
    subroutine test_constant_precision_is_closed_form()
        class(order_prior_offline), allocatable :: sp
        type(image)              :: even, odd, w, ref
        type(pcg_solver_outcome) :: outcome
        real(c_float), pointer   :: rp(:,:,:)
        real, allocatable, target :: x(:,:,:)
        real :: err
        write(logfhandle,'(A)') '-- test_constant_precision_is_closed_form'
        call make_pair(even, odd, .false.)
        allocate(sp)
        call sp%new(even, odd, TMSKDIAM)
        sp%dgrid      = 1.
        sp%data_scale = 1.
        sp%hdiag      = 1.
        call protein_weight(w)
        call ref%new(sp%ldim, sp%smpd)
        call sp%closed_form(1, w, TLAMBDA, ref)
        call sp%even%get_rmat_ptr(rp)
        allocate(x(TBOX,TBOX,TBOX), source=rp(:TBOX,:TBOX,:TBOX))
        call sp%solve_half(1, w, TLAMBDA, x, outcome)
        call ref%get_rmat_ptr(rp)
        err = norm2(x - rp(:TBOX,:TBOX,:TBOX)) / norm2(rp(:TBOX,:TBOX,:TBOX))
        nullify(rp)
        write(logfhandle,'(A,ES10.3,A,I0)') 'relative error ', err, ', CG iterations ', outcome%iteration_count
        call assert_true(outcome%converged, 'the surrogate CG converges with a constant precision')
        call assert_true(err <= (1. + TLAMBDA) * CG_RTOL, &
            &'with a constant precision the surrogate solve is the closed-form shrink h/(h + lambda_s (1-w))')
        call sp%kill
        call even%kill
        call odd%kill
        call w%kill
        call ref%kill
    end subroutine test_constant_precision_is_closed_form

    !> CG needs <A u, v> = <u, A v> and <u, A u> > 0, likewise for the preconditioner
    subroutine test_operator_symmetry()
        class(order_prior_offline), allocatable :: sp
        type(image)        :: even, odd, w
        type(surrogate_op) :: op
        real, allocatable, target :: u(:), v(:), au(:), av(:)
        real(dp) :: auv, uav, uau
        integer  :: nvox
        write(logfhandle,'(A)') '-- test_operator_symmetry'
        call make_pair(even, odd, .true.)
        allocate(sp)
        call sp%new(even, odd, TMSKDIAM)
        call protein_weight(w)
        call sp%init_op(op, w, TLAMBDA)
        nvox = op%size()
        allocate(u(nvox), v(nvox), au(nvox), av(nvox))
        call random_number(u)
        call random_number(v)
        u = u - 0.5
        v = v - 0.5
        call op%apply(u, au)
        call op%apply(v, av)
        auv = op%dot(au, v)
        uav = op%dot(u, av)
        uau = op%dot(u, au)
        call assert_true(abs(auv - uav) <= SYM_TOL * sqrt(op%dot(au,au) * op%dot(v,v)), &
            &'the surrogate operator with the order ridge is symmetric')
        call assert_true(uau > 0._dp, 'the surrogate operator is positive definite on a random vector')
        call op%precond(u, au)
        call op%precond(v, av)
        auv = op%dot(au, v)
        uav = op%dot(u, av)
        uau = op%dot(u, au)
        call assert_true(abs(auv - uav) <= SYM_TOL * sqrt(op%dot(au,au) * op%dot(v,v)), &
            &'the two-region preconditioner is symmetric')
        call assert_true(uau > 0._dp, 'the two-region preconditioner is positive definite on a random vector')
        call sp%kill
        call even%kill
        call odd%kill
        call w%kill
    end subroutine test_operator_symmetry

    !> Noise power rising with k, as in a reconstruction: a smooth micelle lives where D(k) exceeds lambda_s. With an
    !! infinite solvent region the operator passes each shell by D/(D + lambda_s), the shrink by h/(h + lambda_s);
    !! the operator must keep at least half of that predicted gain on the micelle and predict the other half better
    subroutine test_micelle_retention()
        class(order_prior_offline), allocatable :: sp
        type(image)              :: even, odd, w, mic, ce
        type(pcg_solver_outcome) :: outcome
        real(c_float), pointer   :: rp(:,:,:), rm(:,:,:)
        real,     allocatable, target :: xe(:,:,:)
        real,     allocatable    :: spec(:)
        logical,  allocatable    :: lmic(:,:,:)
        real(dp) :: num, den
        real     :: lam_s, ret_cf, ret_op, pred_cf, pred_op
        integer  :: sh
        write(logfhandle,'(A)') '-- test_micelle_retention'
        call make_pair(even, odd, .true., mic)
        allocate(sp)
        call sp%new(even, odd, TMSKDIAM)
        call protein_weight(w)
        lam_s = TLAMBDA * sp%hdiag
        ! the prediction from the construction: the micelle's power per shell weighted by the shell gain
        call mic%fft
        allocate(spec(fdim(TBOX)-1))
        call mic%power_spectrum(spec)
        call mic%ifft
        num = 0._dp
        den = 0._dp
        do sh = 1, size(spec)
            num = num + real(spec(sh),dp) * real(sh,dp)**2 * real(sp%dshell(sh) / (sp%dshell(sh) + lam_s),dp)
            den = den + real(spec(sh),dp) * real(sh,dp)**2
        end do
        pred_op = real(num / den)
        pred_cf = sp%hdiag / (sp%hdiag + lam_s)
        lmic = shell_mask(MIC_R - MIC_W, MIC_R + MIC_W)
        call ce%new(sp%ldim, sp%smpd)
        call sp%closed_form(1, w, TLAMBDA, ce)
        ret_cf = retention(ce)
        call sp%even%get_rmat_ptr(rp)
        allocate(xe(TBOX,TBOX,TBOX), source=rp(:TBOX,:TBOX,:TBOX))
        nullify(rp)
        call sp%solve_half(1, w, TLAMBDA, xe, outcome)
        call ce%set_rmat(xe, .false.)
        ret_op = retention(ce)
        write(logfhandle,'(A,4F8.4)') 'micelle retention shrink/operator, predicted shrink/operator ', &
            &ret_cf, ret_op, pred_cf, pred_op
        call assert_true(pred_op > pred_cf, 'the fixture puts the micelle where the data outweigh the ridge')
        call assert_true(ret_op - ret_cf >= 0.5 * (pred_op - pred_cf), &
            &'the operator keeps the low-resolution density of a smooth micelle that the pointwise shrink removes')
        call sp%kill
        call even%kill
        call odd%kill
        call w%kill
        call mic%kill
        call ce%kill

    contains

        !> regression of the solution on the true micelle inside the micelle shell
        real function retention( sol )
            class(image), intent(inout) :: sol
            real(dp) :: sxy, sxx
            call sol%get_rmat_ptr(rp)
            call mic%get_rmat_ptr(rm)
            sxy = sum(real(rp(:TBOX,:TBOX,:TBOX),dp) * real(rm(:TBOX,:TBOX,:TBOX),dp), mask=lmic)
            sxx = sum(real(rm(:TBOX,:TBOX,:TBOX),dp)**2, mask=lmic)
            retention = real(sxy / sxx)
            nullify(rp, rm)
        end function retention

    end subroutine test_micelle_retention

    !> the threaded loops of the operator and the preconditioner give one thread's result
    subroutine test_team_matches_serial()
        class(order_prior_offline), allocatable :: sp
        type(image)        :: even, odd, w
        type(surrogate_op) :: op
        real, allocatable, target :: u(:), a1(:), a3(:), p1(:), p3(:)
        integer :: nthr_prev, nvox
        write(logfhandle,'(A)') '-- test_team_matches_serial'
        call make_pair(even, odd, .true.)
        allocate(sp)
        call sp%new(even, odd, TMSKDIAM)
        call protein_weight(w)
        call sp%init_op(op, w, TLAMBDA)
        nvox = op%size()
        allocate(u(nvox), a1(nvox), a3(nvox), p1(nvox), p3(nvox))
        call random_number(u)
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call op%apply(u, a1)
        call op%precond(u, p1)
        !$ call omp_set_num_threads(3)
        call op%apply(u, a3)
        call op%precond(u, p3)
        !$ call omp_set_num_threads(nthr_prev)
        call assert_true(maxval(abs(a3 - a1)) <= TEAM_TOL * maxval(abs(a1)), &
            &'a team of three applies the surrogate operator as one thread does')
        call assert_true(maxval(abs(p3 - p1)) <= TEAM_TOL * maxval(abs(p1)), &
            &'a team of three applies the two-region preconditioner as one thread does')
        call sp%kill
        call even%kill
        call odd%kill
        call w%kill
    end subroutine test_team_matches_serial

    !> protein blobs inside PROT_R and a smooth micelle shell, independent noise per half; colored: noise power
    !! rising linearly with k (a reconstruction's), else white
    subroutine make_pair( even, odd, l_colored, mic )
        type(image),           intent(inout) :: even, odd
        logical,               intent(in)    :: l_colored
        type(image), optional, intent(inout) :: mic
        real, parameter :: NOISE_SDEV = 0.6
        real, allocatable :: prot(:,:,:), shell(:,:,:), centres(:,:)
        real    :: r2, d, c
        integer :: i, j, k, b
        call set_fixed_seed(TSEED)
        allocate(centres(3,10))
        call random_number(centres)
        centres = (centres - 0.5) * 2. * (PROT_R - 2.) / sqrt(3.)
        allocate(prot(TBOX,TBOX,TBOX), shell(TBOX,TBOX,TBOX), source=0.)
        c = real(TBOX/2 + 1)
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    do b = 1, size(centres,2)
                        r2 = (real(i) - c - centres(1,b))**2 + (real(j) - c - centres(2,b))**2 + (real(k) - c - centres(3,b))**2
                        prot(i,j,k) = prot(i,j,k) + exp(-0.5 * r2)
                    end do
                    d = sqrt((real(i) - c)**2 + (real(j) - c)**2 + (real(k) - c)**2)
                    shell(i,j,k) = 0.5 * exp(-0.5 * ((d - MIC_R) / MIC_W)**2)
                end do
            end do
        end do
        call even%new([TBOX,TBOX,TBOX], TSMPD)
        call odd%new([TBOX,TBOX,TBOX], TSMPD)
        call even%set_rmat(prot + shell + noise(), .false.)
        call odd%set_rmat(prot + shell + noise(), .false.)
        if( present(mic) )then
            call mic%new([TBOX,TBOX,TBOX], TSMPD)
            call mic%set_rmat(shell, .false.)
        endif

    contains

        function noise() result( rmat )
            real, allocatable :: rmat(:,:,:), filt(:)
            real(c_float), pointer :: rp_noise(:,:,:)
            type(image) :: nimg
            integer     :: sh
            allocate(rmat(TBOX,TBOX,TBOX))
            call random_number(rmat)
            rmat = sqrt(12.) * (rmat - 0.5)
            if( l_colored )then
                call nimg%new([TBOX,TBOX,TBOX], TSMPD)
                call nimg%set_rmat(rmat, .false.)
                allocate(filt(fdim(TBOX)-1))
                do sh = 1, size(filt)
                    filt(sh) = sqrt(real(sh) / real(size(filt)))
                end do
                call nimg%apply_filter(filt)
                call nimg%get_rmat_ptr(rp_noise)
                rmat = rp_noise(:TBOX,:TBOX,:TBOX)
                nullify(rp_noise)
                rmat = rmat / sqrt(sum(rmat**2) / real(TBOX)**3)
                call nimg%kill
            endif
            rmat = NOISE_SDEV * rmat
        end function noise

    end subroutine make_pair

    !> the true protein region, softened over about a voxel
    subroutine protein_weight( w )
        type(image), intent(inout) :: w
        real, allocatable :: rmat(:,:,:)
        real    :: c, d
        integer :: i, j, k
        allocate(rmat(TBOX,TBOX,TBOX))
        c = real(TBOX/2 + 1)
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    d = sqrt((real(i) - c)**2 + (real(j) - c)**2 + (real(k) - c)**2)
                    rmat(i,j,k) = 1. / (1. + exp((d - PROT_R) / 0.7))
                end do
            end do
        end do
        call w%new([TBOX,TBOX,TBOX], TSMPD)
        call w%set_rmat(rmat, .false.)
    end subroutine protein_weight

    function shell_mask( rlo, rhi ) result( l )
        real, intent(in)     :: rlo, rhi
        logical, allocatable :: l(:,:,:)
        real    :: c, d
        integer :: i, j, k
        allocate(l(TBOX,TBOX,TBOX))
        c = real(TBOX/2 + 1)
        do k = 1, TBOX
            do j = 1, TBOX
                do i = 1, TBOX
                    d = sqrt((real(i) - c)**2 + (real(j) - c)**2 + (real(k) - c)**2)
                    l(i,j,k) = d >= rlo .and. d <= rhi
                end do
            end do
        end do
    end function shell_mask

end submodule simple_order_prior_offline_tester
