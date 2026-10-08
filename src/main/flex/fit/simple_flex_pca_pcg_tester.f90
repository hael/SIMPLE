!@descr: white-box self-test of the flex_pca PCG M-step operator (submodule of simple_flex_pca_pcg)
! test_flex_pcg_operator reads flex_pcg_t's private components and kernels, so it is a submodule of the
! operator's module; as a *_tester file it is built only with BUILD_TESTS=ON. The dense-lattice scatters
! and folds it compares the band-list path against (check D) live here with it. simple_flex_pcg_tester
! asserts the six checks by name.
submodule (simple_flex_pca_pcg) simple_flex_pca_pcg_tester
use simple_core_module_api, only: KBWINSZ
use simple_rnd,             only: gasdev, ran3, seed_rnd_fixed
use simple_ori_utils,       only: euler2m
implicit none
#include "simple_local_flags.inc"

contains

    !> Self-test on a random volume: (A) kernel operator vs the exact nonuniform-DFT Gram; (B) 2x rhs deposit vs the
    !! exact adjoint; (C) box <= 32: CG recovery from central-slice samples on a sphere; (D) band-list kernels and
    !! rhs vs the dense fold, bitwise; (E) operator symmetry; (F) unregularised PSD and ridge-regularised PD.
    module subroutine test_flex_pcg_operator( box, nsamples, l_pass, passes, sweep )
        integer,           intent(in)  :: box, nsamples
        logical,           intent(out) :: l_pass
        logical, optional, intent(out) :: passes(6)
        logical, optional, intent(in)  :: sweep       !< (C) over all twelve settings (default: the baseline)
        type(flex_pcg_t) :: op
        type(flex_pcg_outcome_t) :: out
        real,    allocatable :: kacc(:,:), kpk(:,:), u(:,:,:,:), hu(:,:,:,:), eu4(:,:,:,:)
        real,    allocatable :: v4(:,:,:,:), hv4(:,:,:,:), hy4(:,:,:,:), h0x4(:,:,:,:), h0y4(:,:,:,:)
        real,    allocatable :: lam_save(:)
        real,    allocatable :: kacc4(:,:,:,:), kpk4(:,:,:,:)
        complex, allocatable :: racc4(:,:,:,:), rpk4(:,:,:,:)
        real,    allocatable :: locs(:,:), wts(:), b4(:,:,:,:), x4(:,:,:,:), rho_t(:,:,:,:), ut(:,:,:)
        complex, allocatable :: racc(:,:), rpk(:,:)
        real(dp), allocatable :: be(:,:,:), eu(:,:,:), bx(:,:,:)
        complex(dp), allocatable :: ex(:), ey(:), ez(:), ysmp(:)
        complex(dp) :: fval
        complex     :: cv(1)
        real(dp)    :: twopi_n, cc, scale, err, na, nb
        real(dp)    :: xhy, hxy, xhx, yhy, xh0, yh0, sym_rel, sym_tol, lift_rel
        real(dp)    :: xnorm, ynorm, h0xnorm, h0ynorm
        real    :: w(3,3,3), loc(3), loc2(3), ctr, sig, dx, dy, dz, rotp(3,3), ang(3)
        integer :: i, j, k, s, sgn, i0(3), wdim, iwinsz, nyq, nyqsq, ns, ip, h, kk, rho_lb(3), rho_ub(3)
        integer :: c, np, isw
        real    :: mfrac, nlev, lamv, lamr(3)
        real(dp) :: yrms
        logical :: pass_a, pass_b, pass_c, pass_d, pass_e, pass_f, l_sweep
        integer :: tt, ijk(3), nmiss, nsw
        l_sweep = .false.
        if( present(sweep) ) l_sweep = sweep
        call op%new(box, 1.0, 1)
        call op%set_band(op%Rnat)
        wdim   = 2*ceiling(KBWINSZ - 0.5) + 1
        iwinsz = ceiling(KBWINSZ - 0.5)
        if( wdim /= 3 ) THROW_HARD('test_flex_pcg_operator assumes the 3-tap KB window')
        nyq = box/2
        c   = box/2 + 1
        twopi_n = 2.0_dp * PI / real(box,dp)
        allocate(locs(3,nsamples), wts(nsamples))
        call seed_rnd_fixed(20260925)   ! reproducible sample positions, volume modulation and slices
        do s = 1, nsamples
            locs(:,s) = (2.0*[ran3(), ran3(), ran3()] - 1.0) * (0.5*real(nyq))
            wts(s)    = 1.0
        end do
        allocate(u(box,box,box,1), hu(box,box,box,1), eu4(box,box,box,1))
        ctr = real(box)/2.0 + 1.0
        sig = 0.12*real(box)
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    dx = real(i)-ctr; dy = real(j)-ctr; dz = real(k)-ctr
                    u(i,j,k,1) = exp(-(dx*dx+dy*dy+dz*dz)/(2.0*sig*sig)) * (1.0 + 0.3*(ran3()-0.5))
                end do
            end do
        end do
        allocate(ex(box), ey(box), ez(box))
        ! ================= (A) operator against the exact Gram =================
        ! band-limited apodized input Pi(E u), as the operator sees it
        eu4(:,:,:,1) = op%env * u(:,:,:,1)
        call op%bandlimit_img(eu4(:,:,:,1), op%nimg)
        allocate(eu(box,box,box), be(box,box,box), source=0.0_dp)
        eu = real(eu4(:,:,:,1),dp)
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc = real(sgn) * locs(:,s)
                call exps(loc)
                call sample_exact(eu, fval)
                fval = fval * real(wts(s),dp)
                call adjoint_add(fval, be)
            end do
        end do
        ! E Pi (exact Gram): band-limit the reference the way the operator band-limits its output
        eu4(:,:,:,1) = real(be)
        call op%bandlimit_img(eu4(:,:,:,1), op%nimg)
        be = real(op%env,dp) * real(eu4(:,:,:,1),dp)
        ! the kernel operator on the same samples and mates
        call op%alloc_accum(kacc)
        allocate(kacc4(op%npairs, op%lims3(1,2)-op%lims3(1,1)+1, op%lims3(2,2)-op%lims3(2,1)+1, &
            &op%lims3(3,2)-op%lims3(3,1)+1), source=0.0)
        rho_lb = [-(iwinsz+1), -nyq-iwinsz-1, -nyq-iwinsz-1]
        rho_ub = [ nyq+iwinsz+1, nyq+iwinsz+1, nyq+iwinsz+1]
        allocate(rho_t(1, rho_ub(1)-rho_lb(1)+1, rho_ub(2)-rho_lb(2)+1, rho_ub(3)-rho_lb(3)+1), source=0.0)
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc2 = real(op%padf) * real(sgn) * locs(:,s)
                i0   = nint(loc2) - iwinsz
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call op%scatter_pairs_wrap(i0, w, [wts(s)], 1.0, kacc)
                    call scatter_pairs_wrap_dense(op, i0, w, [wts(s)], 1.0, kacc4)
                else
                    call op%scatter_pairs_nowrap(i0, w, [wts(s)], 1.0, kacc)
                    call scatter_pairs_nowrap_dense(op, i0, w, [wts(s)], 1.0, kacc4)
                endif
                call deposit_density(real(sgn)*locs(:,s), wts(s), rho_t)
            end do
        end do
        call op%alloc_packed(kpk)
        call op%fold_accum(kacc, kpk)
        ! (D) band-list accumulation and fold against the dense lattice: bitwise
        allocate(kpk4(op%npairs, op%cdim(1), op%cdim(2), op%cdim(3)), source=0.0)
        call fold_accum_dense(op, kacc4, kpk4)
        nmiss = 0
        do tt = 1, op%npk
            ijk = op%pijk(:,tt)
            if( any(kpk4(:,ijk(1),ijk(2),ijk(3)) /= kpk(:,tt)) ) nmiss = nmiss + 1
            kpk4(:,ijk(1),ijk(2),ijk(3)) = 0.0
        end do
        pass_d = nmiss == 0 .and. all(kpk4 == 0.0)
        write(logfhandle,'(A,I0,A,I0,A,I0,A,L1)') '>>> FLEX PCG TEST (D) band-list kernels vs dense fold: slots=', op%npk, &
            &'  mismatching slots=', nmiss, '  dense mass outside the list=', count(kpk4 /= 0.0), '  pass=', pass_d
        deallocate(kpk4)
        call op%finalize(kpk)
        eu4(:,:,:,1) = op%env * u(:,:,:,1)
        call op%apply_operator(eu4, hu)
        hu(:,:,:,1) = op%env * hu(:,:,:,1)
        call op%prep_floor(rho_t, rho_lb)
        ! Test the B iterated by cg_core directly. The envelope factors above belong only to
        ! the exact-DFT comparison; wrapping just one side in E produces a non-symmetric E*B.
        allocate(v4(box,box,box,1), hv4(box,box,box,1), hy4(box,box,box,1))
        allocate(h0x4(box,box,box,1), h0y4(box,box,box,1))
        do k = 1, box
            do j = 1, box
                do i = 1, box
                    ! Formula vectors leave (C)'s fixed-seed slice stream unchanged.
                    tt = i + box*((j-1) + box*(k-1))
                    eu4(i,j,k,1) = sin(0.173*real(tt)) + 0.5*cos(0.071*real(tt))
                    v4(i,j,k,1)  = cos(0.113*real(tt)) - 0.25*sin(0.053*real(tt))
                end do
            end do
        end do
        allocate(lam_save, source=op%lam)
        op%lam = 0.0
        call op%apply_operator(eu4, h0x4)
        call op%apply_operator(v4, h0y4)
        op%lam = lam_save
        call op%apply_operator(eu4, hv4)
        call op%apply_operator(v4, hy4)
        xhy = sum(real(eu4,dp)*real(h0y4,dp))
        hxy = sum(real(h0x4,dp)*real(v4,dp))
        xhx = sum(real(eu4,dp)*real(hv4,dp))
        yhy = sum(real(v4,dp)*real(hy4,dp))
        xh0 = sum(real(eu4,dp)*real(h0x4,dp))
        yh0 = sum(real(v4,dp)*real(h0y4,dp))
        xnorm  = sqrt(sum(real(eu4,dp)**2)); ynorm  = sqrt(sum(real(v4,dp)**2))
        h0xnorm = sqrt(sum(real(h0x4,dp)**2)); h0ynorm = sqrt(sum(real(h0y4,dp)**2))
        sym_rel = abs(xhy - hxy) / max(1.d0, abs(xhy), abs(hxy))
        sym_tol = 32.d0*real(epsilon(1.0),dp)*sqrt(real(box**3,dp))
        lift_rel = max(abs((xhx-xh0) - real(lam_save(1),dp)*xnorm*xnorm) / &
            &max(1.d0, abs(xhx-xh0)), abs((yhy-yh0) - real(lam_save(1),dp)*ynorm*ynorm) / &
            &max(1.d0, abs(yhy-yh0)))
        pass_e = sym_rel <= sym_tol
        pass_f = xh0 >= -sym_tol*xnorm*h0xnorm .and. yh0 >= -sym_tol*ynorm*h0ynorm .and. &
            &lam_save(1) > 0.0 .and. xhx > 0.d0 .and. yhy > 0.d0 .and. lift_rel <= 10.d0*sym_tol
        write(logfhandle,'(A,ES10.3,A,ES10.3,A,L1)') &
            &'>>> FLEX PCG TEST (E) symmetry: relative error=', sym_rel, ' tolerance=', sym_tol, ' pass=', pass_e
        write(logfhandle,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3,A,L1)') &
            &'>>> FLEX PCG TEST (F) quadratic forms: xH0x=', xh0, ' yH0y=', yh0, &
            &'  actual lam=', real(lam_save(1),dp), ' ridge-lift relative error=', lift_rel, ' pass=', pass_f
        deallocate(v4, hv4, hy4, h0x4, h0y4, lam_save, rho_t)
        call compare3(be, real(hu(:,:,:,1),dp), cc, scale, err, na, nb)
        write(logfhandle,'(A,I0,A,I0,A,F8.5,A,F10.5,A,ES10.3,A,ES10.3)') '>>> FLEX PCG TEST (A) operator box=', box, &
            &' samples=', nsamples, '  kernel vs exact Gram: corr=', real(cc), '  LS scale=', real(scale), &
            &'  rel_resid=', real(err), '  |exact|=', real(na)
        pass_a = abs(scale - 1.0_dp) < 0.05_dp .and. err < 0.1_dp
        ! ================= (B) right-hand side: 2x deposit of exact samples vs exact adjoint =================
        allocate(ysmp(nsamples), source=(0.0_dp,0.0_dp))
        allocate(bx(box,box,box), source=0.0_dp)
        eu = real(u(:,:,:,1),dp)
        do s = 1, nsamples
            loc = locs(:,s)
            call exps(loc)
            call sample_exact(eu, ysmp(s))
        end do
        call op%alloc_rhs_accum(racc)
        allocate(racc4(op%ncomp, op%lims3(1,2)-op%lims3(1,1)+1, op%lims3(2,2)-op%lims3(2,1)+1, &
            &op%lims3(3,2)-op%lims3(3,1)+1), source=cmplx(0.,0.))
        do s = 1, nsamples
            do sgn = 1, -1, -2
                loc  = real(sgn) * locs(:,s)
                loc2 = real(op%padf) * loc
                i0   = nint(loc2) - iwinsz
                if( sgn == 1 )then
                    fval = ysmp(s)
                else
                    fval = conjg(ysmp(s))
                endif
                cv(1) = cmplx(real(real(fval)*real(wts(s),dp), sp), real(aimag(fval)*real(wts(s),dp), sp))
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call op%scatter_rhs_wrap(i0, w, cv, racc)
                    call scatter_rhs_wrap_dense(op, i0, w, cv, racc4)
                else
                    call op%scatter_rhs_nowrap(i0, w, cv, racc)
                    call scatter_rhs_nowrap_dense(op, i0, w, cv, racc4)
                endif
                ! exact adjoint of the same sample
                call exps(loc)
                call adjoint_add(fval*real(wts(s),dp), bx)
            end do
        end do
        call op%alloc_rhs_packed(rpk)
        call op%fold_rhs(racc, rpk)
        allocate(rpk4(op%ncomp, op%cdim(1), op%cdim(2), op%cdim(3)), source=cmplx(0.,0.))
        call fold_rhs_dense(op, racc4, rpk4)
        nmiss = 0
        do tt = 1, op%npk
            ijk = op%pijk(:,tt)
            if( any(rpk4(:,ijk(1),ijk(2),ijk(3)) /= rpk(:,tt)) ) nmiss = nmiss + 1
            rpk4(:,ijk(1),ijk(2),ijk(3)) = cmplx(0.,0.)
        end do
        pass_d = pass_d .and. nmiss == 0 .and. all(rpk4 == cmplx(0.,0.))
        write(logfhandle,'(A,I0,A,I0,A,L1)') '>>> FLEX PCG TEST (D) band-list rhs vs dense fold: mismatching slots=', nmiss, &
            &'  dense mass outside the list=', count(rpk4 /= cmplx(0.,0.)), '  pass=', pass_d
        deallocate(rpk4)
        allocate(b4(box,box,box,1))
        call op%finalize_rhs(rpk, b4)
        eu4(:,:,:,1) = real(bx)
        call op%bandlimit_img(eu4(:,:,:,1), op%nimg)
        bx = real(eu4(:,:,:,1),dp)
        call compare3(bx, real(b4(:,:,:,1),dp), cc, scale, err, na, nb)
        write(logfhandle,'(A,F8.5,A,F10.5,A,ES10.3,A,ES10.3)') '>>> FLEX PCG TEST (B) rhs deposit vs exact adjoint: corr=', &
            &real(cc), '  LS scale=', real(scale), '  rel_resid=', real(err), '  |exact|=', real(na)
        pass_b = abs(scale - 1.0_dp) < 0.05_dp .and. err < 0.1_dp
        deallocate(kpk, rpk, b4, be, bx, ysmp)
        ! ================= (C) preconditioned CG solve from central-slice samples =================
        ! the clean, generously supported baseline (mask 0.40, no noise, no Tikhonov term); the sweep adds
        ! the tighter support, sample noise and the Tikhonov term. Every clean solve must meet the baseline
        ! criterion; the noisy ones characterise what the Tikhonov term is for (unregularised, CG runs
        ! into the null space and the error outside the low-pass grows without bound)
        pass_c = .true.
        if( box <= 32 )then
            np    = 48
            nyqsq = nyq*nyq
            ns = 0
            do ip = 1, np
                do kk = -nyq, nyq
                    do h = -nyq, nyq
                        if( h*h + kk*kk <= nyqsq ) ns = ns + 1
                    end do
                end do
            end do
            deallocate(locs, wts)
            allocate(locs(3,ns), wts(ns))
            wts = 1.0
            s = 0
            do ip = 1, np
                ang(1) = 2.0*PI*ran3(); ang(2) = acos(2.0*ran3()-1.0); ang(3) = 2.0*PI*ran3()
                rotp = euler2m(ang*(180.0/PI))
                do kk = -nyq, nyq
                    do h = -nyq, nyq
                        if( h*h + kk*kk > nyqsq ) cycle
                        s = s + 1
                        locs(:,s) = matmul(real([h,kk,0]), rotp)
                    end do
                end do
            end do
            ! exact samples of the volume
            allocate(ysmp(ns), source=(0.0_dp,0.0_dp))
            eu = real(u(:,:,:,1),dp)
            do s = 1, ns
                call exps(locs(:,s))
                call sample_exact(eu, ysmp(s))
            end do
            yrms = sqrt(sum(abs(ysmp)**2) / real(ns,dp))
            ! kernels on the 2x lattice and the 1x gridding density (independent of the data)
            allocate(rho_t(1, rho_ub(1)-rho_lb(1)+1, rho_ub(2)-rho_lb(2)+1, rho_ub(3)-rho_lb(3)+1), source=0.0)
            call op%alloc_accum(kacc)
            do s = 1, ns
                loc  = locs(:,s)
                loc2 = real(op%padf) * loc
                i0   = nint(loc2) - iwinsz
                call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                if( op%win_wraps(i0) )then
                    call op%scatter_pairs_wrap(i0, w, [1.0], 1.0, kacc)
                else
                    call op%scatter_pairs_nowrap(i0, w, [1.0], 1.0, kacc)
                endif
                call deposit_density(loc, 1.0, rho_t)
            end do
            call op%alloc_packed(kpk)
            call op%fold_accum(kacc, kpk)
            call op%finalize(kpk)
            allocate(b4(box,box,box,1), x4(box,box,box,1), source=0.0)
            allocate(ut(box,box,box), source=0.0)
            nsw = merge(12, 1, l_sweep)
            do isw = 1, nsw
                mfrac = merge(0.4, 0.25, mod((isw-1)/6, 2) == 0)
                nlev  = merge(0.0, 0.5,  mod((isw-1)/3, 2) == 0)
                lamr  = [0.0, 1.0e-3, 1.0e-2]
                lamv  = lamr(mod(isw-1, 3) + 1)
                call op%set_window_sphere(mfrac*real(box))
                call op%set_lambda_relative(lamv)
                ! right-hand sides of the (noisy) samples
                call op%alloc_rhs_accum(racc)
                do s = 1, ns
                    loc  = locs(:,s)
                    loc2 = real(op%padf) * loc
                    i0   = nint(loc2) - iwinsz
                    call op%kbwin%apod_mat_3d_fast(loc2, iwinsz, wdim, w)
                    fval = ysmp(s)
                    if( nlev > 0.0 ) fval = fval + real(nlev,dp)*yrms*cmplx(gasdev(), gasdev(), kind=dp)/sqrt(2.0_dp)
                    cv(1) = cmplx(real(real(fval), sp), real(aimag(fval), sp))
                    if( op%win_wraps(i0) )then
                        call op%scatter_rhs_wrap(i0, w, cv, racc)
                    else
                        call op%scatter_rhs_nowrap(i0, w, cv, racc)
                    endif
                end do
                call op%alloc_rhs_packed(rpk)
                call op%fold_rhs(racc, rpk)
                call op%finalize_rhs(rpk, b4)
                x4 = 0.0
                call op%cg_core(b4, x4, rho_t, rho_lb, 60, 1.0e-4, out, 'FLEX PCG TEST (C)')
                ut = u(:,:,:,1) * op%mask
                err = sqrt(sum((real(x4(:,:,:,1),dp) - real(ut,dp))**2)) / sqrt(sum(real(ut,dp)**2))
                call lowpass(x4(:,:,:,1), nyq/2)
                call lowpass(ut, nyq/2)
                cc = sqrt(sum((real(x4(:,:,:,1),dp) - real(ut,dp))**2)) / sqrt(sum(real(ut,dp)**2))
                write(logfhandle,'(A,F5.2,A,F4.2,A,ES8.1,A,I3,A,I3,A,ES9.2,A,ES9.2,A,ES9.2)') &
                    &'>>> FLEX PCG TEST (C) sweep: mask=', mfrac, ' noise=', nlev, ' lam=', lamv, &
                    &'  iters=', out%iteration_count, '  its_to_1e-2=', out%iters_to_1e2, &
                    &'  final resid=', out%final_rel_residual, '  err full=', real(err), '  err inner=', real(cc)
                if( nlev == 0.0 ) pass_c = pass_c .and. cc < 0.05_dp .and. out%final_rel_residual < 1.0e-2
                deallocate(rpk)
            end do
            write(logfhandle,'(A,I0,A,I0)') '>>> FLEX PCG TEST (C) slices=', np, ' samples=', ns
            deallocate(kpk, b4, x4, rho_t, ut, ysmp)
        endif
        l_pass = pass_a .and. pass_b .and. pass_c .and. pass_d .and. pass_e .and. pass_f
        if( present(passes) ) passes = [pass_a, pass_b, pass_c, pass_d, pass_e, pass_f]
        if( l_pass )then
            write(logfhandle,'(A)') '    PASS: operator, rhs, solve, band lists, symmetry and positivity'
        else
            write(logfhandle,'(A,6L2)') '    FAIL (operator, rhs, solve, band lists, symmetry, positivity): ', &
                &pass_a, pass_b, pass_c, pass_d, pass_e, pass_f
        endif
        call op%kill
        deallocate(locs, wts, u, hu, eu4, eu, ex, ey, ez)

    contains

        !> Add one sample to the native half-lattice density used by prep_floor.
        subroutine deposit_density( p, weight, rho )
            real, intent(in)    :: p(3), weight
            real, intent(inout) :: rho(:,:,:,:)
            real    :: ppos(3), wrho(3,3,3)
            integer :: base(3), ii, jj, ll, ir, jr, lr
            ppos = p
            if( ppos(1) < 0.0 ) ppos = -ppos
            base = nint(ppos) - iwinsz
            call op%kbwin%apod_mat_3d_fast(ppos, iwinsz, wdim, wrho)
            do ll = 1, wdim
                lr = base(3) + ll - rho_lb(3)
                do jj = 1, wdim
                    jr = base(2) + jj - rho_lb(2)
                    do ii = 1, wdim
                        ir = base(1) + ii - rho_lb(1)
                        rho(1,ir,jr,lr) = rho(1,ir,jr,lr) + weight*wrho(ii,jj,ll)
                    end do
                end do
            end do
        end subroutine deposit_density

        !> separable exponentials of one sample position
        subroutine exps( p )
            real, intent(in) :: p(3)
            integer :: ii
            do ii = 1, box
                ex(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(1),dp)*real(ii-c,dp), dp))
                ey(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(2),dp)*real(ii-c,dp), dp))
                ez(ii) = exp(cmplx(0.0_dp, -twopi_n*real(p(3),dp)*real(ii-c,dp), dp))
            end do
        end subroutine exps

        !> F = (1/N^3) sum_n v(n) e^{-i...} for the current exponentials, one axis at a time: the x sum
        !! for every (y,z) column is one real (2,N) x (N,N^2) product (library matmul, also at -O0), then
        !! y and z; (C) takes 38256 central-slice samples at box 32, 1.25e9 voxel terms. v is the box^3
        !! volume by sequence association
        subroutine sample_exact( v, f )
            real(dp),    intent(in)  :: v(box,box*box)
            complex(dp), intent(out) :: f
            real(dp) :: exri(2,box), t(2,box*box)
            integer  :: kk2, j0
            exri(1,:) = real(ex, dp)
            exri(2,:) = aimag(ex)
            t = matmul(exri, v)
            f = (0.0_dp, 0.0_dp)
            do kk2 = 1, box
                j0 = (kk2 - 1)*box
                f = f + ez(kk2) * sum(ey * cmplx(t(1,j0+1:j0+box), t(2,j0+1:j0+box), dp))
            end do
            f = f / real(box,dp)**3
        end subroutine sample_exact

        !> acc(n) += Re[ f e^{+i...} ] for the current exponentials
        subroutine adjoint_add( f, acc )
            complex(dp), intent(in)    :: f
            real(dp),    intent(inout) :: acc(:,:,:)
            integer :: ii, jj, kk2
            !$omp parallel do default(shared) private(ii,jj,kk2) schedule(static)
            do kk2 = 1, box
                do jj = 1, box
                    do ii = 1, box
                        acc(ii,jj,kk2) = acc(ii,jj,kk2) + real(f * conjg(ex(ii)*ey(jj)*ez(kk2)), dp)
                    end do
                end do
            end do
            !$omp end parallel do
        end subroutine adjoint_add

        !> spherical low-pass of a native-lattice volume to shell rmax
        subroutine lowpass( v, rmax )
            real,    intent(inout) :: v(:,:,:)
            integer, intent(in)    :: rmax
            real, pointer :: rp(:,:,:)
            integer :: lims(3,2), hh, kh, mh, ph(3)
            call op%nimg%set_rmat(v, .false.)
            call op%nimg%fft()
            lims = op%nimg%loop_lims(2)
            do mh = lims(3,1), lims(3,2)
                do kh = lims(2,1), lims(2,2)
                    do hh = lims(1,1), lims(1,2)
                        if( nint(sqrt(real(hh*hh + kh*kh + mh*mh))) > rmax )then
                            ph = op%nimg%comp_addr_phys(hh,kh,mh)
                            call op%nimg%set_cmat_at(ph(1),ph(2),ph(3), cmplx(0.,0.))
                        endif
                    end do
                end do
            end do
            call op%nimg%ifft()
            call op%nimg%get_rmat_ptr(rp)
            v = rp(1:box,1:box,1:box)
        end subroutine lowpass

        !> correlation, least-squares scale of b onto a, relative residual, norms
        subroutine compare3( a, bb, cc_, scale_, err_, na_, nb_ )
            real(dp), intent(in)  :: a(:,:,:), bb(:,:,:)
            real(dp), intent(out) :: cc_, scale_, err_, na_, nb_
            real(dp) :: num, den
            na_ = sqrt(sum(a**2))
            nb_ = sqrt(sum(bb**2))
            num = sum(a*bb)
            den = sum(bb**2)
            cc_ = num / max(na_*nb_, 1.0e-30_dp)
            scale_ = 1.0_dp
            if( den > 0.0_dp ) scale_ = num/den
            err_ = sqrt(sum((a - scale_*bb)**2)) / max(na_, 1.0e-30_dp)
        end subroutine compare3

    end subroutine test_flex_pcg_operator

    ! ---------------- dense-lattice references of check (D) ----------------

    pure subroutine scatter_pairs_nowrap_dense( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        real    :: wv
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    wv = w(di,dj,dk) * val
                    kacc(:,ih,ik,im) = kacc(:,ih,ik,im) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_nowrap_dense

    pure subroutine scatter_pairs_wrap_dense( self, i0, w, dpack, val, kacc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim), dpack(:), val
        real,              intent(inout) :: kacc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        real    :: wv
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    wv = w(di,dj,dk) * val
                    kacc(:,ih,ik,im) = kacc(:,ih,ik,im) + wv * dpack(:)
                end do
            end do
        end do
    end subroutine scatter_pairs_wrap_dense

    pure subroutine scatter_rhs_nowrap_dense( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        do dk = 1, self%wdim
            im = i0(3) + dk - 1 - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = i0(2) + dj - 1 - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = i0(1) + di - 1 - self%lims3(1,1) + 1
                    racc(:,ih,ik,im) = racc(:,ih,ik,im) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_nowrap_dense

    pure subroutine scatter_rhs_wrap_dense( self, i0, w, vals, racc )
        class(flex_pcg_t), intent(in)    :: self
        integer,           intent(in)    :: i0(3)
        real,              intent(in)    :: w(self%wdim,self%wdim,self%wdim)
        complex,           intent(in)    :: vals(:)
        complex,           intent(inout) :: racc(:,:,:,:)
        integer :: di, dj, dk, ih, ik, im
        do dk = 1, self%wdim
            im = self%wrap(i0(3)+dk-1) - self%lims3(3,1) + 1
            do dj = 1, self%wdim
                ik = self%wrap(i0(2)+dj-1) - self%lims3(2,1) + 1
                do di = 1, self%wdim
                    ih = self%wrap(i0(1)+di-1) - self%lims3(1,1) + 1
                    racc(:,ih,ik,im) = racc(:,ih,ik,im) + w(di,dj,dk) * vals(:)
                end do
            end do
        end do
    end subroutine scatter_rhs_wrap_dense

    subroutine fold_accum_dense( self, kacc, kpk )
        type(flex_pcg_t),  intent(inout) :: self
        real, allocatable, intent(inout) :: kacc(:,:,:,:)
        real,              intent(inout) :: kpk(:,:,:,:)
        integer :: h, hh, k, m, phys(3), ih, ik, im
        if( .not. allocated(kacc) ) return
        if( size(kpk,1) /= self%npairs ) THROW_HARD('packed kernel set has the wrong leading extent; fold_accum')
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,phys,ih,ik,im) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, -self%lims3(1,1)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    ih   = hh - self%lims3(1,1) + 1
                    ik   = k  - self%lims3(2,1) + 1
                    im   = m  - self%lims3(3,1) + 1
                    kpk(:,phys(1),phys(2),phys(3)) = kpk(:,phys(1),phys(2),phys(3)) + kacc(:,ih,ik,im)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(kacc)
    end subroutine fold_accum_dense

    subroutine fold_rhs_dense( self, racc, rpk )
        type(flex_pcg_t),     intent(inout) :: self
        complex, allocatable, intent(inout) :: racc(:,:,:,:)
        complex,              intent(inout) :: rpk(:,:,:,:)
        integer :: h, hh, k, m, phys(3), ih, ik, im
        if( .not. allocated(racc) ) return
        if( size(rpk,1) /= self%ncomp ) THROW_HARD('packed rhs set has the wrong leading extent; fold_rhs')
        !$omp parallel do collapse(2) default(shared) private(h,hh,k,m,phys,ih,ik,im) schedule(static)
        do m = self%lims3(3,1), self%lims3(3,2)
            do k = self%lims3(2,1), self%lims3(2,2)
                do h = 0, -self%lims3(1,1)
                    hh   = self%wrap(h)
                    phys = self%wimg%comp_addr_phys(h,k,m)
                    ih   = hh - self%lims3(1,1) + 1
                    ik   = k  - self%lims3(2,1) + 1
                    im   = m  - self%lims3(3,1) + 1
                    rpk(:,phys(1),phys(2),phys(3)) = rpk(:,phys(1),phys(2),phys(3)) + racc(:,ih,ik,im)
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(racc)
    end subroutine fold_rhs_dense

end submodule simple_flex_pca_pcg_tester
